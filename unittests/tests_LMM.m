rng(1)
input = [];
for k = 1:30
     input{k} = simulate_data_lmm_v2('noise',0.1,'noise_components',1,'srate',20,'n_events',64)

end
%%

EEG = um_designmat(input,'inputGroupingName','subject',...
    'eventtypes','sim',...
    'formula','y~1+cat(condA)*cat(condB)+(1+condA*condB|subject)+(1|stimulus)','codingschema','effects');

EEG= um_timeexpandDesignmat(EEG,'timelimits',[-0.1,0.5]);
%%
% Currently I recommend the bobyqa optimizer. Seems to be faster
model_fv = um_mmfit(EEG,input,'channel',63,'optimizer','fminunc','covariance','CompSymm'); % Todo: Directly read covariance from formula like lme4

%%

model_fv = model_fv.initstats;
results = fixedEffects(model_fv,0.05,'satterthwaite');

%%
  model = model_fv
    R = model.Psi.NumBlocks;
    covtable = cell(R,1);
    covmat = covtable;
    offset = 0;
    tbl = covarianceParameters(model,0.05,0);
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
 
        offset = endk;
    end
    colorbar
    
    figure,
    nTimeshifts = model.Psi.NumBlocks;
    subplot(2,nTimeshifts,1:nTimeshifts);
    plot(reshape(model.betaHat,nTimeshifts,[]),'o-')
    for k = 1:R
        subplot(2,nTimeshifts,nTimeshifts + k);
        imagesc(covmat{k})
        caxis(prctile(tbl.Estimate,[1,99]))
    end