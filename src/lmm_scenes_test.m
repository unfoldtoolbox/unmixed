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

EEG = um_designmat(input,'eventtypes','fixation','formula','y~1+spl(sac_amp_incoming,10)+(1+spl(sac_amp_incoming,10)|subject)');
EEG= um_timeexpandDesignmat(EEG,'timelimits',[-0.1,0.5]);

%%
model_fv = um_mmfit(EEG,input,'channel',1,'optimizer','bobyqa','covariance','Diagonal');

save(['model_scene_' DataHash(model_fv)],'model_fv')
% model_fminunc, model_q
% 
% model = model_fminunc
%% Plotting
if 1 == 0
    %%
    model = model_fv
    R = model.Psi.NumBlocks;
    covtable = cell(R,1);
    covmat = covtable;
    offset = 0;
    tbl = covarianceParameters(model,0.05,0)
    for k = 1:R
        % (1) Build an index idxk to extract relevant rows from tbl
        % for grouping variable k.
        matk = model.Psi.Matrices{k};
        startk = offset + 1;
        endk = offset + matk.NumParametersExcludingSigma;
        idxk = startk : endk;
        
        % (2) k the element of covtable.
        covtable{k} = [tbl(idxk,:)];
        covmat{k}= triu(ones(ceil(sqrt(length(idxk)))));
        covmat{k}(covmat{k}==1) = tbl.Estimate(idxk);
        % (3) Add a title for covtable{k}.
        %ttl = ['Group: ',model.GroupingInfo.GNames{k},...
        %    ', Covariance Type: ',model.slme.Psi.Matrices{k}.Type];
        % <entry key="Title_covtable">Covariance Type: {0}</entry>
        %     ttl = getString(message('stats:LinearMixedModel:Title_covtable',model.slme.Psi.Matrices{k}.Type));
        %     covtable{k} = classreg.regr.lmeutils.titleddataset(covtable{k},ttl);
        
        % (4) Update offset to go to the next grouping variable.
        offset = endk;
    end
    colorbar
    
    figure,
    nTimeshifts = model.Psi.NumBlocks;
    subplot(2,nTimeshifts,1:nTimeshifts);
    plot(reshape(model.betaHat,nTimeshifts,2),'o-')
    for k = 1:R
        subplot(2,nTimeshifts,nTimeshifts + k);
        imagesc(covmat{k})
        caxis(prctile(tbl.Estimate,[1,99]))
    end
end