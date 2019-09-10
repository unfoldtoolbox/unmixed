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


result = cfg;
result.events = length(EEG.event);
result.sizeXdc = size(EEG.unmixed.uf_fixef.Xdc);
result.sizeZdc = size(EEG.unmixed.uf_ranef{1}.Zdc);
tic
result.model = um_mmfit(EEG,input,'channel',cfg.channel,'optimizer',cfg.optimizer);
result.timing = toc;
toc
save(fullfile(cfg.folder,DataHash(cfg)),'result','cfg','input')
