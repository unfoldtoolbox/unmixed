function [EEG] = simulate_data_lmm_v2(varargin)

simCFG= finputcheck(varargin,...
    {'n_events','integer',[],100; 
    'epochlength','real',[],0.5; %in s
    'simulationtype','string',{'realistic','ideal','ideal_hanning'},'realistic'; %realistic uses SEREEGA, ideal generates a single channel respons
    
    'noise_components','real',[],10; % number of random noise components with brown noise
    'noise','real',[],1; % strength of noise
    'u_noise','real',[],0.5; %variability of noise level - subjectwise
    
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
    'overlapparam','',[],[[0.35,0.1,0,0];[0.1,0,0,0]]; % Effect Coding with -1 / 1! [interceptmean, diffCondA,diffCondB; interceptSD diffCondA ...] for lognormal
                                                           %['interceptMin,interceptMax] ...
    'overlapminimum','real',[],0.1; %whats the shortest time two stimuli can follow eachother? This will be adjusted AFTER taking the previous parameters into acount, thus biasing them!
    
    },'mode','error');

assert(~ischar(simCFG),simCFG)




%%  Generate data
sig.time = simCFG.srate * simCFG.epochlength;

switch simCFG.simulationtype
    case 'realistic'
        simulated_data = um_sereega_epochs('n_epochs',simCFG.n_events,...
            'noise_components',simCFG.noise_components,...
            'srate',simCFG.srate,'noise_orient',1,'epochlength',simCFG.epochlength);
    case 'ideal'
        % warning currently p1 p3 and n1 all occur at the same time
        simulated_data.p1.data(1,:,:) = repmat([1 zeros(1,sig.time-1)]',1,simCFG.n_events);
        simulated_data.p3.data(1,:,:) = repmat([1 zeros(1,sig.time-1)]',1,simCFG.n_events);
        simulated_data.n1.data(1,:,:) = repmat([1 zeros(1,sig.time-1)]',1,simCFG.n_events);
        simulated_data.random.data(1,:,:) = randn(1,size(simulated_data.n1.data,2),size(simulated_data.n1.data,3)*3);
    case 'ideal_hanning'
        % warning currently p1 p3 and n1 all occur at the same time
        simulated_data.p1.data(1,:,:) = repmat(hanning(sig.time),1,simCFG.n_events);
        simulated_data.p3.data(1,:,:) = zeros(sig.time,simCFG.n_events);%repmat([1 zeros(1,sig.time-1)]',1,simCFG.n_events);
        simulated_data.n1.data(1,:,:) = zeros(sig.time,simCFG.n_events);%repmat([1 zeros(1,sig.time-1)]',1,simCFG.n_events);
        simulated_data.random.data(1,:,:) =  randn(1,size(simulated_data.n1.data,2),size(simulated_data.n1.data,3)*3);
end
%%Generate Stimulus Timings
% How much time should the continuous EEG have
%   ~ 1 stim / s
% whatTimeForEvents = simCFG.srate*simCFG.epochlength*simCFG.n_events*2;
% howManyEvents = ceil(simCFG.n_events*1.5);
% 
% while true
%     
%     switch simCFG.overlaptype
%         case 'uniform'
%             % we should also use the cumsum approach (see below)
% %             warning('you might need to adapt whatTimeForEvents, or wait quite long ;-)')
% %             ix = unique(sort(randi(whatTimeForEvents,howManyEvents,1)));
%             ix = cumsum(round(simCFG.srate*(simCFG.overlapparam(1) + rand(howManyEvents,1)*simCFG.overlapparam(2))));
%         case 'lognormal'
%             
%             m =  simCFG.overlapparam(1)*simCFG.srate;% on average events should start at ~250ms
%             v = (simCFG.overlapparam(2)*simCFG.srate)^2; % we want ~100ms jitter, gives nice distributions
%             % taken from 'lognstat' matlab help:
%             mu = log((m^2)/sqrt(v+m^2));
%             sigma = sqrt(log(v/(m^2)+1));
%             % to visualize:
%             % figure,hist(lognrnd(mu,sigma,1,1000)/EEG.srate,1000)
%             ix = cumsum(round(lognrnd(mu,sigma,1,howManyEvents)));
%     end
%     ix = round(ix);
%     del = diff(ix) <simCFG.overlapminimum*simCFG.srate; % min overlap should be 100ms
%     ix(del) = [];
%     if length(ix) >= howManyEvents/1.5
%         ix = ix(1:howManyEvents/1.5);
%         break
%     end
% end
% ix = ix+3*simCFG.srate*simCFG.epochlength; % add a bit of slack in the beginning :-)
% assert(length(ix)==length(EEG.event))

%% Generate Designmatrix
% Subject =  1*rand()   + 0.5*slope_rand() +  2*rand() + 2  *slope_rand()

% Intercept, Main Effect A, Main Effect B, Interaction
X = [ones(1,simCFG.n_events);%  intercept
    2*randi([0,1],1,simCFG.n_events)-1;%  factor e.g. trialtype
    2*randi([0,1],1,simCFG.n_events)-1]';
X(:,end+1) = X(:,2) .* X(:,3); % interaktion


adjustedOverlap = X*simCFG.overlapparam';
% generate overlap
ix = nan(1,simCFG.n_events);
for ev = 1:simCFG.n_events
    ix(ev) = -1; % to keep the while running
   while ix(ev) < simCFG.overlapminimum*simCFG.srate
        switch simCFG.overlaptype
        case 'uniform'
            
            ix(ev) = round(simCFG.srate*(adjustedOverlap(ev,1) + rand(1)*(adjustedOverlap(ev,2)-adjustedOverlap(ev,1))));
        case 'lognormal'
            
            m =  adjustedOverlap(ev,1)*simCFG.srate;% on average events should start at ~250ms
            v = (adjustedOverlap(ev,2)*simCFG.srate)^2; % we want ~100ms jitter, gives nice distributions
            % taken from 'lognstat' matlab help:
            mu = log((m^2)/sqrt(v+m^2));
            sigma = sqrt(log(v/(m^2)+1));
            % to visualize:
            % figure,hist(lognrnd(mu,sigma,1,1000)/EEG.srate,1000)
            ix(ev) = cumsum(round(lognrnd(mu,sigma,1,1)));
        end
        
    end
    


end
ix = round(cumsum(ix));


stimA_0 = mod(1:sum(X(:,2)==-1),simCFG.n_items)+1;
stimA_1 = mod(1:sum(X(:,2)==1 ),simCFG.n_items)+1+max(stimA_0);
stimA_0 = stimA_0(randperm(length(stimA_0)));
stimA_1 = stimA_1(randperm(length(stimA_1)));

stimA = nan(length(X(:,2)),1);
stimA(X(:,2)==-1) = stimA_0;
stimA(X(:,2)==1)  = stimA_1;

EEG = eeg_emptyset();
EEG.event = struct('type','sim',...
    'latency',num2cell(ix'),...
    'trialnum',num2cell(1:simCFG.n_events)',...
    'condA',num2cell(X(:,2)/2+0.5),...
    'condB',num2cell(X(:,3)/2+0.5),...
    'stimulus',num2cell(stimA))';
EEG.data = [];
EEG.epoch = [];
EEG.srate = simCFG.srate;
EEG.trials = 1;

% EEG=uf_designmat(EEG,'formula','y~1+condA+condB','eventtypes','sim');
% EEG.pnts =ceil(EEG.event(end).latency + simCFG.epochlength*simCFG.srate); % fix (prospective) timing for last event

EEG.data = zeros(size(simulated_data.p1.data,1),ceil(max(ix)+1.2*simCFG.srate*simCFG.epochlength));
% EEG2 = uf_timeexpandDesignmat(EEG,'timelimits',[0,simCFG.epochlength]);
for fn = fieldnames(simulated_data)'
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
%         warning('same item variability for all 3 peaks')
        rng(old);
        item_effect = item_effect_nitems(stimA); % what each item has in addtion
%         X(:,1) = 1 + item_effect/(b(1)+Z(1)); % I simply add this to each intercept
        
        
    end
    
  
    simdat_raw = simulated_data.(fn{1}).data;
    simdat = nan(size(simdat_raw));
    for ch = 1:size(simdat,1)
        multWith = (X*(b+Z)')';
        if simCFG.randomItem
            multWith = multWith + item_effect';
        end
    simdat(ch,:,:) =bsxfun(@times,squeeze(simdat_raw(ch,:,:)),multWith);

    end
    
    for k = 1:length(EEG.event)
        EEG.data(:,ix(k):(ix(k)+simCFG.srate*simCFG.epochlength-1)) = EEG.data(:,ix(k):(ix(k)+simCFG.srate*simCFG.epochlength-1))+simdat(:,:,k);
    
    end
    EEG.sim.X = X;
    EEG.sim.Z.(fn{1}) = Z;
    EEG.sim.u.(fn{1}) = u;
    EEG.sim.b.(fn{1}) = b;
    EEG.sim.uCov = D*R*D;
    if simCFG.randomItem
        EEG.sim.item = item_effect;
    end
end
EEG.sim.simCFG = simCFG;
%% adding noise
% simCFG.noise = 0.01;
% To abs or not to abs: I think it actually does not matter (because sign of noise is arbitrary), but I like it
% if the noise parameter is positive
if simCFG.noise >0
noiselevel = abs(randn(1)*simCFG.u_noise+simCFG.noise);
if size(EEG.data,2) > size(simulated_data.random.data(1,:),2)
    
   warning('Overlap is not large enough, need to duplicating noise') 
   % in case the "olverapped" EEG is larger than n_events*epochlength, we
   % do not have enough epochs to span the whole continuous overlapped EEG
   % and thus need to copy. Times two should be enough for most cases
   % TODO: Make a principled noise function, or generate more noise
   % epochs, or regenerate noise epochs?!
   simulated_data.random.data(1,end:size(simulated_data.random.data(1,:),2)) = simulated_data.random.data(:,:);
end
noise_norm = simulated_data.random.data(:,1:size(EEG.data,2));
scale_noise_by = prctile(noise_norm(:),[10,90]);
noise_norm = noise_norm./diff(scale_noise_by);
EEG.data = EEG.data + noiselevel*noise_norm;

EEG.sim.noiselevel = noiselevel;
else 
    EEG.sim.noiselevel = 0;
end
if isfield(simulated_data.p1,'chanlocs')
    EEG.chanlocs = simulated_data.p1.chanlocs;
end
EEG = eeg_checkset(EEG,'eventconsistency');

if 1 == 0
    %% For Debugging purposes, plot the results
    EEG2=uf_designmat(EEG,'formula','y~1+cat(condA)*cat(condB)','eventtypes','sim','codingschema','effects');
    EEG2 = uf_timeexpandDesignmat(EEG2,'timelimits',[0,simCFG.epochlength]);
    EEG2 = uf_glmfit(EEG2);
    
    uf_plotParam(uf_condense(EEG2),'channel',63)
end