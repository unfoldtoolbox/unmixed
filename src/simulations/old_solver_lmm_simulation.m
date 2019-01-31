
cfg = struct();
cfg.itemeffect = 0;
X = [];

Y = [];
G1 = [];
G2 = [];
G1_real = {};
G2_real = {};
beta_dc = [];
beta_nodc = [];
rng(1)
for k = 1:52
    
    EEG = simulate_data_lmm('noise',5,'srate',20,'basis','dirac');
    EEG.pnts = size(EEG.data,2);
    
    EEG = uf_designmat(EEG,'eventtypes','sim','formula','y~1+b');
    %     EEG = uf_designmat(EEG,'eventtypes','sim','formula','y~1+slope+trialtype+contrast');
    timelim = [0,0.25];
    EEG = uf_timeexpandDesignmat(EEG,'timelimits',timelim);
%     EEG = uf_glmfit(EEG);
%     whereToAdd = EEG.unfold.Xdc(:,EEG.unfold.Xdc_terms2cols==1); % nsamp x nlocalsamp
    if cfg.itemeffect
        
        
        itemeffect= cumsum(whereToAdd);
        
        randomcomponent = [0;rand(max(itemeffect(:)),1)*500];
        
        
        
        itemeffect(whereToAdd==0) = 0;
        %     itemeffect_full = itemeffect + randomcomponent(1+(full(itemeffect)));
        itemeffect_full = 0 + randomcomponent(1+(full(itemeffect)));
        
        
        whatToAdd = [EEG.unfold.times>=0 &EEG.unfold.times<1];
        %     EEG.data = EEG.data + (itemeffect_full * whatToAdd')';
        
        %     G2 = cat(1,G2,full(cumsum(whereToAdd(:,EEG.unfold.times == 0)))); %overlap possible!
        G2_real{end+1} = itemeffect'; % staircase matrix
        
    end
    % add trial specific effect
    X = cat(1,X,EEG.unfold.Xdc);
    Y = cat(1,Y,EEG.data(1,:)');
    %     G1 = cat(1,G1,ones(size(EEG.unfold.Xdc,1),1)*k); %no overlap possible
    G1_real{end+1} = EEG.unfold.Xdc;
    
%     EEG = uf_epoch(EEG,'timelimits',timelim);
%     EEG = uf_glmfit_nodc(EEG);
%     beta_dc = cat(2,beta_dc, EEG.unfold.beta_dc(:));
%     beta_nodc = cat(2,beta_nodc, EEG.unfold.beta_nodc(:));
    
end



%%

% nTimeshifts    = size(G1_real{1},1);
% G1_real_concat = cat(2,G1_real{:})';
% if cfg.itemeffect
%     G2_real_concat = cat(2,G2_real{:})';
%     nItems = length(G2_real);
% end
% Zall =

% Zs = sparse(size(G1_real_concat,1), size(G1_real_concat,2) .* nSubjects);
% for k = 1:nSubjects
%     if k == 0
%         continue
%     end
%     ix = ((k-1)*size(G1_real{1},1)+1):((k)*size(G1_real{1},1));
%     Zs(:,ix) = Zs(:,ix) + G1_real_concat==k;
% end
% 
% if cfg.itemeffect
%     Zi = sparse(size(G1_real_concat,1), size(G1_real_concat,2) .* nItems);
%     
%     for k = 1:nItems
%         if k == 0
%             continue
%         end
%         ix = ((k-1)*size(G1_real{1},1)+1):((k)*size(G1_real{1},1));
%         Zi(:,ix) = Zi(:,ix) + G2_real_concat==k;
%     end
    
    
% end

% Bring timeshifted Zsubject matrix in correct form (not 100% sure this is
% the correct way, but somehow mtalab needs to know what effects belong
% together)
% 
%%
nTimeshifts = length(EEG.unfold.times);
nSubjects =length(G1_real);

% All Random Slopes
Zs = repmat(X,1,nSubjects);

% Zs = Zs >0;
P = 1:size(Zs,2);
P = reshape(P,nTimeshifts,[],nSubjects);

P = permute(P,[2 3 1]);


% P2 = 1:size(P,2);
% P2 = reshape(P2,[],nSubjects);
% P2 = P2';
% P = P(:,P2);
% P = P';
Zs_perm = Zs(:,P(:));
if cfg.itemeffect
    Z = cat(2,Zs_perm,Zi);
else
    Z = Zs_perm;
end



FitMethod = 'ML';


% Generate Psi
% gnames = {'g1','g2'};
% ZColNames = {{'G11'},{'G11'}};

if cfg.itemeffect
    Glevels = [repmat(nSubjects,1,nTimeshifts),repmat(nSubjects,1,nTimeshifts)];
else
    Glevels = [repmat(nSubjects,1,nTimeshifts)];
end
mat = [];
for i = 1:length(Glevels)
    mat{i} = classreg.regr.lmeutils.covmats.CovarianceMatrix.createCovariance('Full',2);
    %                    1,...model.RandomInfo.q(i),
    %                    'Name',gnames{i},... model.GroupingInfo.GNames{i},...
    %                     'VariableNames',ZColNames{i});%model.RandomInfo.ZColNames{i});
end

Psi = classreg.regr.lmeutils.covmats.BlockedCovariance(mat,Glevels);

% Psi = classreg.regr.lmeutils.covmats.DiagonalCovariance(1);% is actually: classreg.regr.lmeutils.covmats.BlockedCovariance
dofit = true;
dostats = false;


%%
ix = any(X,2);
tic
addpath('src\lmm\temporaryFunctions\')
model = classreg.regr.lmeutils.StandardLinearMixedModel(X(ix,:),Y(ix),Z(ix,:),Psi,FitMethod,dofit,dostats,'OptimizerOptions',struct('Display','Iter'));
% b = classreg.regr.lmeutils.StandardLinearMixedModel(X,Y,Z,Psi,FitMethod,dofit,dostats,'Optimizer','fminunc','OptimizerOptions',optimoptions('fminunc','Display','Iter'))%struct('Display','Iter'));
rmpath('src\lmm\temporaryFunctions\')
toc
%%
R = length(Glevels);
covtable = cell(R+1,1);
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

%% Plot
figure,

subplot(2,nTimeshifts,1:nTimeshifts);
plot(reshape(model.betaHat,nTimeshifts,2),'o-')
for k = 1:R
    subplot(2,nTimeshifts,nTimeshifts + k);
    imagesc(covmat{k})
    caxis(prctile(tbl.Estimate,[1,99]))
end


