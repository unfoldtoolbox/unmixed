function lmm_scenes_test()
if 1 == 0
    cd '/home/predatt/benehi/projects/unmixed'
    init_unmixed
    if 1 == 0
        %%
        addpath('/home/common/matlab/fieldtrip/qsub')
        qsubfeval(@lmm_scenes_test,'memreq',10*1024^3,'timreq',60*47*60)
    end
end
%%
input = [];
d  = dir(fullfile('raw','scenes','*.set'));
input = cellfun(@(x,y)fullfile(x,y),{d.folder},{d.name},'uniformoutput',0);

for k = 1:length(input)
    EEG = scenes_1st_preprocessing(k,d(k).folder);
    EEG =pop_select(EEG,'channel',42);
    EEG = pop_resample(EEG,50);
    input{k}=EEG;
end
%%

EEG = um_designmat(input,'eventtypes','fixation','formula','y~1+spl(sac_amp_incoming,10)+(1|subject)');
EEG= um_timeexpandDesignmat(EEG,'timelimits',[-0.1,0.5]);

%%
model_fv = um_mmfit(EEG,input,'channel',1,'optimizer','bobyqa','covariance','FullCholesky');

save(['model_scene_' DataHash(model_fv)],'model_fv')
