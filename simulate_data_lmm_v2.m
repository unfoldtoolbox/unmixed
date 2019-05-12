function [EEG] = simulate_data_lmm_v2(varargin)

simCFG= finputcheck(varargin,...
    {'n_events','integer',[],100; 
    'epochlength','real',[],0.5; %in s
    
    'noise_components','real',[],10; % number of random noise components with brown noise
    'noise','real',[],1; % strength of noise
    
    'b_p1_2x2','real',[],[10,5,-1,3]; % P1: Intercept, MainA, MainB, Inter - effect coded beta
    'u_p1_2x2','real',[],[5,2,2,2]; %   P1: Subject variability
    
    'b_p3_2x2','real',[],[6,3,0,0]; 
    'u_p3_2x2','real',[],[3,2,1,1];
    
    'b_n1_2x2','real',[],[-6,0,-2,-2];
    'u_n1_2x2','real',[],[3,1,1,1];
    
    'srate','integer',[],100; % sampling rate
    'randomItem','boolean',[],0; %add an item effect?
    'u_p1_item','real',[],[5]; % Item effect strength
    'u_p3_item','real',[],[1];
    'u_n1_item','real',[],[4];
    'n_items','real',[],10;... % how many different items per factor 
    'overlaptype','string',{'uniform','lognormal'},'lognormal'; % overlap between events
    },'mode','ignore');

assert(~ischar(simCFG),simCFG)




%% Generate Stimulus Timings
% How much time should the continuous EEG have
%   ~ 1 stim / s
sereega_data = um_sereega_epochs('n_epochs',simCFG.n_events,...
    'noise_components',simCFG.noise_components,...
    'srate',simCFG.srate,'noise_orient',1,'epochlength',simCFG.epochlength);
%%
whatTimeForEvents =simCFG.srate*simCFG.epochlength*simCFG.n_events*0.8;
howManyEvents = simCFG.n_events;

while true
    
    switch simCFG.overlaptype
        case 'uniform'
            % we should also use the cumsum approach (see below)
            warning('you might need to adapt whatTimeForEvents, or wait quite long ;-)')
            ix = unique(sort(randi(whatTimeForEvents,howManyEvents,1)));
            
        case 'lognormal'
            
            m = simCFG.srate*0.25;% on average events should start at ~250ms
            v = (0.1*simCFG.srate)^2; % we want ~100ms jitter, gives nice distributions
            % taken from 'lognstat' matlab help:
            mu = log((m^2)/sqrt(v+m^2));
            sigma = sqrt(log(v/(m^2)+1));
            % to visualize:
            % figure,hist(lognrnd(mu,sigma,1,1000)/EEG.srate,1000)
            ix = cumsum(round(lognrnd(mu,sigma,1,howManyEvents)));
    end
    ix = round(ix);
    del = diff(ix) <0.1*simCFG.srate; % min overlap should be 100ms
    ix(del) = [];
    if length(ix) == howManyEvents
        break
    end
end
ix = ix+3*simCFG.srate*simCFG.epochlength; % add a bit of slack in the beginning :-)
% assert(length(ix)==length(EEG.event))

%% Generate Designmatrix
% Subject =  1*rand()   + 0.5*slope_rand() +  2*rand() + 2  *slope_rand()

% Intercept, Main Effect A, Main Effect B, Interaction
X = [ones(1,length(ix));%  intercept
    2*randi([0,1],1,length(ix))-1;%  factor e.g. trialtype
    2*randi([0,1],1,length(ix))-1]';
X(:,end+1) = X(:,2) .* X(:,3); % interaktion

stimA_0 = mod(1:sum(X(:,2)==-1),simCFG.n_items)+1;
stimA_1 = mod(1:sum(X(:,2)==1 ),simCFG.n_items)+1+max(stimA_0);
stimA_0 = stimA_0(randperm(length(stimA_0)));
stimA_1 = stimA_1(randperm(length(stimA_1)));

stimA = nan(length(X(:,2)),1);
stimA(X(:,2)==-1) = stimA_0;
stimA(X(:,2)==1)  = stimA_1;

EEG = eeg_emptyset();
EEG.event = struct('type','sim',...
    'latency',num2cell(ix)',...
    'trialnum',num2cell(1:length(ix))',...
    'condA',num2cell(X(:,2)/2+0.5),...
    'condB',num2cell(X(:,3)/2+0.5),...
    'stimulus',num2cell(stimA))';
EEG.data = [];
EEG.epoch = [];
EEG.srate = simCFG.srate;
EEG.trials = 1;

% EEG=uf_designmat(EEG,'formula','y~1+condA+condB','eventtypes','sim');
% EEG.pnts =ceil(EEG.event(end).latency + simCFG.epochlength*simCFG.srate); % fix (prospective) timing for last event

EEG.data = zeros(size(sereega_data.p1.data,1),whatTimeForEvents);
% EEG2 = uf_timeexpandDesignmat(EEG,'timelimits',[0,simCFG.epochlength]);
for fn = fieldnames(sereega_data)'
    if strcmp(fn{1},'random')
        continue
    end
    eval(sprintf('b = simCFG.b_%s_2x2;',fn{1}));
    eval(sprintf('u = simCFG.u_%s_2x2;',fn{1}));
    % this adds variability per subject
    % Main Effect A and intercept correlated, Main Effect B and Interaction
    % Correlated :shrug:
    R = [1   0.5   0   0;
        0.5 1     0   0;
        0   0     1   0.5;
        0   0     0.5 1];
    D = diag(u);
    Z = mvnrnd([0,0,0,0],D*R*D);
    
    if simCFG.randomItem
        eval(sprintf('u_item = simCFG.u_%s_item;',fn{1}));
        
        old = rng(1); % to get same effect of each stimulus for each subject in consecutive simulation-runs
        item_effect_nitems = randn(simCFG.n_items*2,1)*u_item;
        warning('same item variability for all 3 peaks')
        rng(old);
        item_effect = item_effect_nitems(stimA); % what each item has in addtion
        X(:,1) = 1 + item_effect/(b(1)+Z(1)); % I simply add this to each intercept
        
        
    end
    
  
    simdat_raw = sereega_data.(fn{1}).data;
    simdat = nan(size(simdat_raw));
    for ch = 1:size(simdat,1)
    simdat(ch,:,:) =bsxfun(@times,squeeze(simdat_raw(ch,:,:)),(X*(b+Z)')');
    end
    
    for k = 1:length(EEG.event)
        EEG.data(:,ix(k):(ix(k)+simCFG.srate*simCFG.epochlength-1)) = EEG.data(:,ix(k):(ix(k)+simCFG.srate*simCFG.epochlength-1))+simdat(:,:,k);
    
    end
    EEG.sim.X = X;
    EEG.sim.Z.(fn{1}) = Z;
    EEG.sim.u.(fn{1}) = u;
    EEG.sim.b.(fn{1}) = b;
    if simCFG.randomItem
        EEG.sim.item = item_effect;
    end
end
%%
simCFG.noise = 0.01;
EEG.data = EEG.data + simCFG.noise*sereega_data.random.data(:,1:size(EEG.data,2));

EEG.chanlocs = sereega_data.p1.chanlocs;
EEG = eeg_checkset(EEG,'eventconsistency');

if 1 == 0
    %% For Debugging purposes, plot the results
EEG2=uf_designmat(EEG,'formula','y~1+cat(condA)*cat(condB)','eventtypes','sim','codingschema','effects');
EEG2 = uf_timeexpandDesignmat(EEG2,'timelimits',[0,simCFG.epochlength]);
EEG2 = uf_glmfit(EEG2);

uf_plotParam(uf_condense(EEG2),'channel',63)
end