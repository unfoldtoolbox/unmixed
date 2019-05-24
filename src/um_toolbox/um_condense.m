function umresult = um_condense(EEG,varargin)
cfg = finputcheck(varargin,...
    {'pvaluesFixed','boolean',[],0;... # calculate parametric pvalues using Satterhwaite DFs
    'pvaluesRandom','boolean',[],0;... # calculate parametric pvalues using Satterhwaite DFs
    },'mode','error');


model = EEG.unmixed.modelfit;
umresult = struct();

%% Fixed Effects
if cfg.pvaluesFixed
    fprintf('Calculating Fixed Effects DF with Satterhwaite approximation')
    umresult.fixef= fixedEffects(model,0.05,'satterthwaite');
else
    fixef = fixedEffects(model,0.05,'residual');
    fixef(:,{'pValue','DF','Lower','Upper'}) = array2table(nan(size(fixef,1),4));
end
nEffects = length(EEG.unmixed.uf_fixef.colnames);
nTimes = length(EEG.unmixed.uf_fixef.times);

fixefmatrix= reshape(table2array(fixef),nEffects,nTimes,7);
umresult.fixef = table();
umresult.fixef.names = EEG.unmixed.uf_fixef.colnames';
for n = 1:length(fixef.Properties.VariableNames)
    name = fixef.Properties.VariableNames{n};
    umresult.fixef.(lower(name))= fixefmatrix(:,:,n);
end
umresult.times = EEG.unmixed.uf_fixef.times;

%% Random Effects


if cfg.pvaluesRandom
    %     fprintf('Calculating Random Effects DF with Satterhwaite approximation')
    %     raneflist = randomEffects(model_fv,0.05,'satterthwaite');
    error('check if this is equal order to the one below')
    % also compare to
    raneflist = covarianceParameters(model,0.05,1); % consists of variances, covariances & last entry is error
else
    raneflist = covarianceParameters(model,0.05,0); % consists of variances, covariances & last entry is error
end

%% Extract the Random Covariance Matrices
% the names are build like this "g:NAME-t:0.01", looking for the dash to
% split

nameSplit = cellfun(@(x)strsplit(x.Name,'-t:'),model.Psi.Matrices,'UniformOutput',0);
nameSplit = vertcat(nameSplit{:});
randomtimes = cellfun(@(x)str2num(x),nameSplit(:,2));
randomgroups = cellfun(@(x)x(3:end),nameSplit(:,1),'UniformOutput',0);

groupingFactors = unique(randomgroups,'stable');

counterIdx = 0;
for groupIdx = 1:length(groupingFactors)
    % (1) Build an index idxk to extract relevant rows from tbl
    % for grouping variable k.
    
    ranef = [];
    % Which rows belong to that groupIdx?
    grIdx = strcmp(groupingFactors{groupIdx},randomgroups);
    
    % We keep track of the time so we can build the time x N x N
    % covmatrices
    timeIdx = 1;
    for rIdx = find(grIdx)'
        matk = model.Psi.Matrices{rIdx};
        startk = counterIdx + 1;
        endk = counterIdx + matk.NumParametersExcludingSigma;
        idxk = startk : endk;
        
        % (2) k the element of covtable.
        %     covtable{rIdx} = [raneflist(idxk,:)];
        
        triangleMatrix = triu(ones(ceil(sqrt(length(idxk)))));

        %init cov
        if strcmp(model.Psi.Matrices{rIdx}.Type,'Diagonal')
            covmat = eye(size(triangleMatrix));
            covmat(covmat==1) = raneflist.Estimate(idxk);
            ranef.corrmat(timeIdx,:,:) = eye(size(covmat));
        else
            covmat= triangleMatrix;
            %fill the estimates
            covmat(triangleMatrix==1) = raneflist.Estimate(idxk);
            % fill the estimates, lower triangle
            covmat(triangleMatrix'==1) = raneflist.Estimate(idxk);
            
            if any(isnan(covmat(:)))
                warning('Nans found in ',model.Psi.Matrices{rIdx}.Name, 'model mightnot be converged')
                ranef.corrmat(timeIdx,:,:) = nan(size(covmat));
            else
                try
                    ranef.corrmat(timeIdx,:,:) = corrcov(covmat);
                catch
                    warning('Covariance Matrix is not positive semidefinite, model might not be converged')
                    ranef.corrmat(timeIdx,:,:) = nan(size(covmat));
                end
            end    
            
        
            
        end
        ranef.covmat(timeIdx,:,:)  = covmat;
        counterIdx = endk;
        timeIdx = timeIdx + 1;
        ranef.groupingvariableLevels=  model.Psi.NumReps(rIdx);
        ranef.variableNames = model.Psi.Matrices{rIdx}.VariableNames;
        
        ranef.type = model.Psi.Matrices{rIdx}.Type;
        ranef.name = groupingFactors{groupIdx};
        
        umresult.ranef(groupIdx) = ranef;
    end
    umresult.logLikelihood = EEG.unmixed.modelfit.loglikHat;
    
end
