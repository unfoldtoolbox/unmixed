classdef um_CustomLinearMixedModel  < classreg.regr.lmeutils.StandardLinearMixedModel
    methods (Access=protected)
        
        function [thetaHat,cause] = doMinimization(slme,fun,theta0)
            
            iscorr = cellfun(@(x)strcmp(x.Type,'corr'),slme.Psi.getCanonicalParameterNames,'UniformOutput',0);
            opts = slme.OptimizerOptions;
            opts.isCovariance = cat(1,iscorr{:})';
            [thetaHat] = ...
                um_bobyqa(fun,theta0,opts);
            
            cause = 0; % XXX Find Exitflags of bobyqa to give more reasonable causes
        end
    end
end