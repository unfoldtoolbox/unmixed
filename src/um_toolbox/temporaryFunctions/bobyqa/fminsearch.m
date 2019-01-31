function  [thetaHat,tmp,cause] = fminsearch(fun,theta0,varargin)
fprintf('Hacked fminqn to bobyqa function! \n') 

%need to do this because of weird functionname
% tmpfun = @(x)fun(x);

% profile on
% fun(theta0)
% profile off
[thetaHat, ~] = bobyqa(fun, theta0,0,1e20,struct('display', 'iter','maxFunEval',50000));

cause = 1;
tmp = [];tmp2 = [];