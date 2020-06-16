function um_simulations_grid(varargin)
%% 
% This script runs simulations on the donders-grid. The starting script 
[cfg] = finputcheck(varargin,...
    {'timelimits','real',[],[-.2 1.2]
    'noise','real',[],20
    'u_noise','real',[],20
    'srate','real',[],50
    'n_events','real',[],10
    'n_subjects','real',[],50
    'optimizer','string',[],'bobyqa'
    'randomItem','','',1;
    'simulationtype','','','ideal_hanning';
    'overlaptype','','','lognormal';
    'b_p1_2x2','','',[5 5 0 0];
    'u_p1_2x2','','',[10 1 0 0];
    'b_p3_2x2','','',[0 0 0 0];
    'u_p3_2x2','','',[0 0 0 0];
    'b_n1_2x2','','',[0 0 0 0];
    'u_n1_2x2','','',[0 0 0 0];
    'u_p3_item','','',0;
    'u_n1_item','','',0;
    'u_p1_item','','',0;
    'noise_components','','',5;
    'folder','','','cache';
    'channel','','',1;
    'formula','','','y~1+condA+condB+(1+condA|subject)';
    'covariance','string',[],'diagonal'
    },'mode','ignore');

% cfg.folder = 'cache_BobyqaAndQausinewton';

if ~exist(cfg.folder,'dir')
    mkdir(cfg.folder);
end
% dont calculate again if already calculated
if exist(fullfile('cache',[DataHash(cfg) '.mat']),'file')
    return
end


rng(1)
input = [];
for k = 1:cfg.n_subjects    
    input{k} = simulate_data_lmm_v2(cfg);
end


EEG = um_designmat(input,'eventtypes','sim','formula',cfg.formula);
EEG= um_timeexpandDesignmat(EEG,'timelimits',cfg.timelimits);


result_lmm = cfg;
result_lmm.events = length(EEG.event);
result_lmm.sizeXdc = size(EEG.unmixed.uf_fixef.Xdc);
result_lmm.sizeZdc = size(EEG.unmixed.uf_ranef{1}.Zdc);
tic
result_lmm.EEG = um_mmfit(EEG,input,'channel',cfg.channel,'optimizer',cfg.optimizer,'covariance',cfg.covariance);
result_lmm.timing = toc;
toc

% 2-stage, I should probably write a function for this at some point
% :shrug:

result_twostage = [];
for subj = 1:length(input)
    EEGsing = uf_designmat(input{subj},'formula',char(EEG.unmixed.formulaFixef),'eventtypes','sim',...
        'codingschema','reference');
    EEGsing = uf_timeexpandDesignmat(EEGsing,'timelimits',cfg.timelimits);
    EEGsing = uf_glmfit(EEGsing,'channel',cfg.channel);
    % 2 stage without overlapcorrection
    EEGsing = uf_epoch(EEGsing,'timelimits',cfg.timelimits);
    EEGsing = uf_glmfit_nodc(EEGsing);
    ufresult = uf_condense(EEGsing);
    if isempty(result_twostage)
        result_twostage = ufresult;
        result_twostage.subject = subj;
    else
        result_twostage.unfold(end+1) = ufresult.unfold;
        result_twostage.beta(:,:,:,end+1) = ufresult.beta;
        result_twostage.beta_nodc(:,:,:,end+1) = ufresult.beta_nodc;
        result_twostage.subject(end+1) = subj;
    end
end

% complete pooling
result_completePooling = fitlm(result_lmm.EEG.unmixed.modelfit.X,result_lmm.EEG.unmixed.modelfit.y);

save(fullfile(cfg.folder,DataHash(cfg)),'result_lmm','cfg','input','result_twostage','result_completePooling')
