function  [thetaHat,tmp,cause] = fminsearch(fun,theta0,varargin)
fprintf('Hacked fminqn to bobyqa function! \n') 

%need to do this because of weird functionname
% tmpfun = @(x)fun(x);
[thetaHat, ~] = bobyqa(fun, theta0,[],[],struct('display', 'iter','maxFunEval',100000));

cause = 1;
tmp = [];tmp2 = [];