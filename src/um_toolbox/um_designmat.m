function [EEG,EEG_fixef] = um_designmat(input,varargin)
[cfg, ufdesignmatCFG] = finputcheck(varargin,...
    {'formula','string',[],[];...
    'inputGroupingName','string',[],'subject';...
    },'mode','ignore');

if ischar(cfg)
error(cfg)
end
%%
cfg.formula = regexprep(cfg.formula,'[\s]','');
ranefRegexp = '\+\(([-\*a-zA-z\(\)0-9\,+]+)\|([a-zA-Z\(\)\,0-9]+)\)';


cfg.formulaRanef= regexp(cfg.formula,ranefRegexp,'tokens');

for k = 1:length(cfg.formulaRanef)
    % add y~ in front for uf_designmat
    cfg.formulaRanef{k}{1} = strcat('y~',cfg.formulaRanef{k}{1});
end
tmp = regexp(cfg.formula,ranefRegexp,'split');
cfg.formulaFixef = sprintf("%s",tmp{:});

% Make sure a grouping variable is not used as an effect
 assert(any(cellfun(@(x)any(strfind(cfg.formulaFixef,x{2})),cfg.formulaRanef))==0,'Please dont add grouping variables also as fixed effects')

%% Load Data info

if iscell(input)
    if isstruct(input{1})
        % nothing to do here
        
    elseif ischar(input)
        % read input files
        for k = 1:length(input)
            input{k} = pop_loadset('filename',input{k},'loadmode','info');
            input{k}.data = []; 
        end
    end
else
    error('either give a cellarray of EEG structures or a cell array of filename.set-strings')
end
%%
%Check equal sampling rate
assert(length(unique(cellfun(@(x)x.srate,input))) == 1,'Not all sets have the same sampling rate')

EEG = eeg_emptyset();
% collect events

max_prev_latency = 0;
% Go through each "subject" and assign a unique subject variable
% Also change the event latencies, so that we can concatenate the designs
for k = 1:length(input)
    singleEEG =  input{k}.event;
    % incase the inputGroupingName already exists, give a warning before
    % overwriting with 1:length().
    if ~(all(cellfun(@(x)isempty(x),strfind(fieldnames(singleEEG),cfg.inputGroupingName))))
        warning('overwritting singleEEG.event.%s with k=1:length(input)',cfg.inputGroupingName)
    end
    
    for e = 1:length(singleEEG)
        
        singleEEG(e).(cfg.inputGroupingName) = k;
        singleEEG(e).urlatency = singleEEG(e).latency;
        singleEEG(e).latency = singleEEG(e).latency + max_prev_latency;
        
    end
    
    EEG.event =[EEG.event singleEEG];
    EEG.pnts = EEG.pnts + input{k}.pnts; % needed for appropriate uf_timeexpand call
    
    max_prev_latency = max_prev_latency + input{k}.pnts;%max([singleEEG(:).latency]);
end




fprintf('um_designmat: Modeling Fixed Effect Part \n')
% there can be only one fixef effects part
EEG_fixef =uf_designmat(EEG,'formula',char(cfg.formulaFixef),ufdesignmatCFG{:});
EEG.unmixed.uf_fixef = EEG_fixef.unfold;

% clear EEG_fixef


fprintf('um_designmat: Modeling Random Effects Part\n')
% but there can be many ranef effect parts
for k = 1:length(cfg.formulaRanef)
    fprintf('um_designmat: Grouping Variable: %s\n',cfg.formulaRanef{k}{2})
    % adding the grouping variables to the designmatrix so that we can
    % timeexpand them later
    formula = [char(cfg.formulaRanef{k}{1}) '+' [cfg.formulaRanef{k}{2}]];
%     formula = char(cfg.formulaRanef{k}{1});
    
    EEG_ranef =uf_designmat(EEG,'formula',formula,ufdesignmatCFG{:});
    EEG_ranef.unfold.ranefgrouping = cfg.formulaRanef{k}{2};
    
    % change variabletype to ranefgrouping
    ix = strcmp(EEG_ranef.unfold.variablenames, cfg.formulaRanef{k}{2});
    EEG_ranef.unfold.variabletypes{ix} = 'ranefgrouping';
    
    if EEG_ranef.unfold.codingschema =="effects"
        % add back the mean
        ix = EEG_ranef.unfold.cols2variablenames(end);
        EEG_ranef.unfold.X(:,ix) =EEG_ranef.unfold.X(:,ix)+EEG_ranef.unfold.effects_mean(end);
        EEG_ranef.unfold.effects_mean(end) = nan;
    end
    % Save it
    EEG.unmixed.uf_ranef{k} = EEG_ranef.unfold;
    EEG.unmixed.datapoints_readin = cellfun(@(x)x.pnts,input); % I might not need this one
    
    
end

clear EEG_ranef
EEG.unmixed.formula = cfg.formula;
EEG.unmixed.formulaFixef = cfg.formulaFixef;
EEG.unmixed.formulaRanef =  cfg.formulaRanef;
EEG.srate = input{k}.srate;

%%
