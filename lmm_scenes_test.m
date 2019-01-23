
input = [];
d  = dir('C:\Users\behinger\scenes\*.set');
input = cellfun(@(x,y)fullfile(x,y),{d.folder},{d.name},'uniformoutput',0);

for k = 1:length(input)
    EEG = scenes_1st_preprocessing(k,'C:\Users\behinger\scenes');
    input{k}= pop_select(EEG,'channel',42);
end
%%

EEG = um_designmat(input,'eventtypes','fixation','formula','y~1+sac_amp_incoming+(1+sac_amp_incoming|subject)');
EEG= um_timeexpandDesignmat(EEG,'timelimits',[-0.1,0.5]);

%%
model_fv = um_mmfit(EEG,input,'channel',1,'optimizer','quasinewton')
% model_fminunc, model_q
model = model_fv
% model = model_fminunc
%% Plotting
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
