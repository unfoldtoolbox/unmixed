function  [thetaHat,tmp,cause] = um_bobyqa(fun,theta0,varargin)
fprintf('Hacked fminqn to bobyqa function! \n') 

%need to do this because of weird functionname
% tmpfun = @(x)fun(x);

% profile on
% fun(theta0)
% profile off
% varargin.MaxFunEvals

%set defaults for variance parameters
ub = repmat(1e20,length(theta0),1);
lb = repmat(-1e20,length(theta0),1);

% warning(' RHO END TO e-2, change back to e-6!')

[thetaHat, fVal] = bobyqa(fun, theta0,lb,ub,struct('display', 'iter',...
    'rho_end',1e-6,'maxFunEval',varargin{1}.MaxFunEvals));

cause = 1;
tmp = [];tmp2 = [];