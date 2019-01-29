function EEG = scenes_1st_preprocessing(subject,varargin)
if nargin == 1
    scenepath = '/net/store/nbp/users/behinger/projects/deconvolution/raw/scenes/';
else
    scenepath = varargin{1};
end
try
    EEG = pop_loadset(fullfile(scenepath,sprintf('MSDARK3_%02i_ica_synched.set',subject)));
    
catch
    error('subject not found')
end

%% preprocessing
EEG = eeg_checkset(EEG);

% rename the separated fixations of the different condition to the same
% name. We will model the effect of condition directly in the formula
for e = 1:length(EEG.event)
    if ismember(EEG.event(e).type,{'fixation1','fixation2','fixation3'})
        EEG.event(e).type = 'fixation';
    end
end

% the angle is not added to fixation, we need to do this manually
if ~exist('allfix','var')
    load(fullfile(scenepath,'fixation_table.mat'));
end

% compied tis from Olaf. first part selects subject, second part only
% condition 1+2+3
ix_sub = allfix(:,22) == subject & ismember(allfix(:,24),[1 2 3]);
subfix = allfix(ix_sub,:);
% first saccade has no incoming saccade
ix_no_incoming_sac = isnan(subfix(:,47));
subfix(ix_no_incoming_sac,:)=[];

% check that EEG saccades and table-saccades fit
assert(size(subfix,1) == sum(strcmp({EEG.event.type},'fixation')),'error, fixation file and EEG file do not match')

% find out which events are fixations
fixIX = find(strcmp({EEG.event.type},'fixation'));


% generate predictors relating
angle = subfix(:,46);
for e = 1:length(fixIX)
    ev = fixIX(e);
    EEG.event(ev).angle = angle(e);
    EEG.event(ev).angle_cos1f = cos(angle(e));
    EEG.event(ev).angle_sin1f = sin(angle(e));
    % observation: whenever we use the second frequency set, the non-unfold
    % estimates explode.
    EEG.event(ev).angle_cos2f = cos(2*angle(e));
    EEG.event(ev).angle_sin2f = sin(2*angle(e));
end

buttonix = find(strcmp({EEG.event.type},'S 99'));
targetix = find(strcmp({EEG.event.type},'S  1'));
stimix = (strcmp({EEG.event.type},'S 12') | strcmp({EEG.event.type},'S 13') | strcmp({EEG.event.type},'S 11'));
darkcondix =  ismember(find(stimix),find(strcmp({EEG.event.type},'S 11')));

for e = 1:length(darkcondix)
    
    EEG.event(targetix(e)).scene_condition = 1;
    if subject == 6 && e>38
        % in this subject, one buttonpress is missing
        EEG.event(targetix(e)).scene_condition = 1;
        e = e-1;
    end
    EEG.event(buttonix(e)).scene_condition = 1;
    
    
end



EEG.event([EEG.event(:).scene_condition] == 1) = []; % delete events that are recorded in darkness
EEG.event([EEG.event(:).scene_condition] == 3) = []; % delete events that are recorded scrambled
% EEG = pop_resample(EEG,250);


EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG);