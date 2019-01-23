function [EEG] = simulate_data_lmm(varargin)

simCFG= finputcheck(varargin,...
    {'datasamples','integer',[],[];
    'noise','real',[],1;
    'srate','integer',[],10;
    'randomItem','boolean',[],0;
    'basis','string',{'box','hanning','dirac','posneg'},'box'
    },'mode','ignore');

assert(~ischar(simCFG),simCFG)

if isempty(simCFG.datasamples)
    simCFG.datasamples = 60*simCFG.srate;
end

sig = struct();
sig.time= 0:1/simCFG.srate:(1-1/simCFG.srate); % 1 second stimulus

sig.shape=zeros(1,length(sig.time),1);
switch simCFG.basis
    case 'box'
        sig.shape= ones(1,length(sig.time));
    case 'hanning'
        sig.shape = hanning(length(sig.time)); %P3
    case 'dirac'
        sig.shape = [0 1 zeros(1,length(sig.time)-2)];
    case 'posneg'
        sig.shape = [hanning(floor(length(sig.time)/2)); -hanning(ceil(length(sig.time)/2))]'; %P3
    otherwise
        error('unknown shape')
end


%% Generate Stimulus Timings
%   ~ 1 stim / s
whatTimeForEvents = simCFG.datasamples-length(sig.shape)*3;
howManyEvents = simCFG.datasamples/simCFG.srate*1;
while true
    ix = unique(sort(randi(whatTimeForEvents,howManyEvents,1)));
%     del = diff(ix) <2;
%     ix(del) = [];
    if length(ix) == howManyEvents
        break
    end
end
ix = ix+5*(1+length(sig.shape)); %

%% Generate Designmatrix
% Subject =  1*rand()   + 0.5*slope_rand() +  2*rand() + 2  *slope_rand()

u = [15,10,0,0];
b = [20,-5,0,0];

% this adds variability per subject
R = [1   0.5   0   0;
    0.5 1     0   0;
    0   0     1   0.5;
    0   0     0.5 1];
D = diag(u);
Z = mvnrnd([0,0,0,0],D*R*D);

%Z = Z.*0;
% this adds main effects
X = [ones(1,length(ix));%  intercept
    rand(1,length(ix)); % slope subject effect, e.g. difficulty
    randi([0,1],1,length(ix));%  factor e.g. trialtype
    rand(1,length(ix));]';% another slope e.g. contrast
warning('WARNING: PUT MOST EFFECTS TO 0')


if simCFG.randomItem
    b(5) = 1;
    Z(5) = 0;
    
   
    old = rng(1); % to get same mean number for itemeffects in every call
    item = randn(size(X,1),1)+5; 
    rng(old);
    X(:,5) = randn(size(X,1),1)*2+item; % some variation around the true mean
    

end



EEG = eeg_emptyset();
EEG.event = struct('type','sim',...
    'latency',num2cell(ix),...
    'trialnum',num2cell(1:length(ix))',...
    'a',num2cell(X(:,1)),...
    'b',num2cell(X(:,2)),...
    'c',num2cell(X(:,3)),...
    'd',num2cell(X(:,4)))';

EEG.srate = simCFG.srate;
if simCFG.randomItem
    EEG=uf_designmat(EEG,'formula','y~-1+a+b+c+d+trialnum','eventtypes','sim');
else
    EEG=uf_designmat(EEG,'formula','y~-1+a+b+c+d','eventtypes','sim');
end




EEG.pnts =EEG.event(end).latency + length(sig.shape);
EEG2 = uf_timeexpandDesignmat(EEG,'timelimits',[0,1]);


tmpB = repmat(sig.shape,length(b),1);
tmpB = bsxfun(@times,tmpB',b+Z);

EEG.data = (EEG2.unfold.Xdc*tmpB(:))';
EEG.data(:,(end+1):(end+length(sig.shape)*10)) = 0;

EEG.data = EEG.data + simCFG.noise * randn(size(EEG.data));
EEG.sim = struct('shape',sig.shape);
EEG = rmfield(EEG,'unfold');
EEG.sim.X = X;
EEG.sim.Z = Z;
EEG.sim.u = u;
EEG.sim.b = b;
EEG.pnts = size(EEG.data,2);
EEG.srate = simCFG.srate;
% EEG = eeg_checkset(EEG,'eventconsistency');