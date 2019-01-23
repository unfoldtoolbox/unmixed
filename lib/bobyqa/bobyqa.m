function [x, fval] = bobyqa(fun, x0, lb, ub, options)
% The BOBYQA algorithm for bound constrained optimization without
% derivatives by M.J.D. Powell
% 
% 
% ==== License ====
% 
% Copyright (c) [2019] [Karlsruhe Institute of Technology
%                       Institute of Engineering Mechanics]
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following condition:
% 
%   * The above copyright notice and this permission notice shall be
%     included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
% 
% ==== Preparation ====
% 
% This file uses the dlib C++ implementation of BOBYQA.
% To use the algorithm, the mex file has to be compiled for the specific
% archtecture of your computer. Please refer to the official MATLAB help
% for further details.
% 
% You need three files to run the algorithm:
% 
%   * This file, which defines default parameters and executes the mex file.
%   * The file "bobyqasub.m" which updates the screen during the iterations
%     and enables the C++ algorithm to evaluate MATLAB objective functions.
%   * The file "bobyqa_alg.cpp" which contains the C++ source code defining
%     a gateway between MATLAB and the dlib library containing the
%     algorithm.
%   * The dlib library containing the C++ source files with the BOBYQA
%     algorithm. Place a copy of the dlib-folder in the same directory as the
%     other three files or adjust the include command for mex function
%     compilation accordingly.
%     You can download the dlib library at http://dlib.net
% 
% The function will display the necessary commands for compilation in case
% of an error.
% 
% 
% ==== Execution ====
% 
% This file sets defaults for the BOBYQA parameters (if not specified) and
% executes the mex function. The mex function itself calls this file to
% evaluate the MATLAB objective function during the iterations.
% 
% Input parameters:
% 
%   * fun           string/function_handle   name of the objective function
% 
%   * x0            vector                   initial optimization vector
% 
%   * lb            vector (optional)        lower bounds for optimization
%                                               default: -1e100 (no bound)
% 
%   * ub            vector (optional)        upper bounds for optimization
%                                               default:  1e100 (no bound)
% 
%   * options       structure                options for display/algorithm
%           .display         display mode:
%                              'none': no output during iterations
%                              'iter': (default) output after every step
%           .npt             number of points for quadratic approximation
%                              default: 2*n + 1
%           .rho_beg         initial trust region radius
%                              default: 10
%           .rho_beg         final trust region radius
%                              default: 1e-6
%           .maxFunEval      maximum number of function evaluations
%                              default: 1000
% 
% Output parameters:
% 
%   * x             vector   final optimization vector
% 
%   * fval          scalar   value of objective function for final step
% 
% If the optional parameters are not specified, the default values are set.
% If only one of the bounds (lb or ub) is specified, the default is set for
% both bounds!
%
%
% ==== Example ====
% 
% To test the algorithm you can run the following MATLAB code as example:
% 
%     fprintf('\n\n>>> First run:\n\n');
%     testfun = @(y) norm(y-[3;5;1;7]);
%     x0 = [-4;5;99;3];
%     lb = -1e100*ones(4,1);
%     ub =  1e100*ones(4,1);
%     options = struct('display', 'iter', 'npt', 9, 'rho_beg', 10, ...
%                      'rho_end', 1e-6, 'maxFunEval', 100);
%     [x, fval] = bobyqa(testfun, x0, lb, ub, options)
%     fprintf('\n\n>>> Second run:\n\n');
%     options = struct('display', 'none', 'npt', 9, 'rho_beg', 10, ...
%                      'rho_end', 1e-6, 'maxFunEval', 1000);
%     [x, fval] = bobyqa(testfun, x0, lb, ub, options)

if (nargin == 0)      % output: compilation info
  dirdlib = dir(pwd);
  dirdlib = dirdlib([dirdlib.isdir]);
  dirdlib = dirdlib(arrayfun(@(info) numel(info.name), dirdlib) > 4);
  if isempty(dirdlib) || strcmpi(dirdlib, 'dlib')
    dirdlib = '';
  else
    dirdlib = dirdlib(find(arrayfun(@(info) strcmpi(info.name(1:4), 'dlib'), dirdlib), 1, 'first')).name;
    dirdlib = [', strcat(''-I"'',pwd,''\\',dirdlib,'"'')'];
  end
  
  dirbob = mfilename('fullpath');
  tmp    = strfind(dirbob,'\');
  dirbob = strrep(extractBefore(dirbob,tmp(end)), '\', '\\');

  fprintf(['\n\n', 'To compile the BOBYQA_ALG mex function, navigate to', '\n\n' ...
           '    ', dirbob, '\n\n', ...
           'and use the command', '\n\n', ...
           '    mex(strcat(''-I"'',pwd,''"'')',dirdlib,', ''bobyqa_alg.cpp'')', '\n\n']);
  
  x = [];
  fval = [];
  
else                  % normal execution
  
  % Number of variables
  n = numel(x0);
  
  % Set options
  opts = struct('display',    'iter', ...      % default options
                'npt',        2*n+1, ...
                'rho_beg',    10,     ...
                'rho_end',    1e-6,   ...
                'maxFunEval', 1000); 
  if (nargin == 5)   % user supplied options
    cOpts = {'display', 'npt', 'rho_beg', 'rho_end', 'maxFunEval'};
    for i = 1:numel(cOpts)
      if isfield(options, cOpts{i})
          opts.(cOpts{i}) = options.(cOpts{i});
      end
    end
  end
  
  % Default bounds
  if ~exist('lb','var') || isempty(lb), lb = -1e100*ones(n,1); end
  if ~exist('ub','var') || isempty(ub), ub =  1e100*ones(n,1); end
  
%   if isa(fun, 'function_handle'), fun = func2str(fun); end
  
  % First call to "bobyqasub.m" => initialization
  bobyqasub(x0, fun, opts);
  
  % Start iteration
  try
    [fval, x] = bobyqa_alg(x0, n, opts.npt, lb, ub, opts.rho_beg, opts.rho_end, opts.maxFunEval);
    bobyqasub(x);       % reset "bobyqasub.m"
  catch ME %#ok<NASGU>
    ME
    dirdlib = dir(pwd);
    dirdlib = dirdlib([dirdlib.isdir]);
    dirdlib = dirdlib(arrayfun(@(info) numel(info.name), dirdlib) > 4);
    if isempty(dirdlib) || strcmpi(dirdlib, 'dlib')
      dirdlib = '';
    else
      dirdlib = dirdlib(find(arrayfun(@(info) strcmpi(info.name(1:4), 'dlib'), dirdlib), 1, 'first')).name;
      dirdlib = [', strcat(''-I"'',pwd,''\\',dirdlib,'"'')'];
    end
    
    dirbob = mfilename('fullpath');
    tmp    = strfind(dirbob,'\');
    dirbob = strrep(extractBefore(dirbob,tmp(end)), '\', '\\');
    
    fprintf(['\n\n', 'An error occured trying to evaluate the BOBYQA_ALG mex function.', '\n', ...
             'If the error perisits try to recompile the mex function for your system.', '\n', ...
             'Navigate to', '\n\n' ...
             '    ', dirbob, '\n\n', ...
             'and use the command', '\n\n', ...
             '    mex(strcat(''-I"'',pwd,''"'')',dirdlib,', ''bobyqa_alg.cpp'')', '\n\n']);
    
    x = [];
    fval = [];
  end
end

end