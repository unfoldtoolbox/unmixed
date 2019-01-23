function fval = bobyqasub(x, fun, options)
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
% ==== Description ====
% 
% This file uses the dlib C++ implementation of BOBYQA.
% To use the algorithm, the mex file has to be compiled for the specific
% archtecture of your computer. Please refer to the official MATLAB help
% for further details.
% 
% The file uses a persistent variable to save previous results for display
% purposes. If a critical error occurs during the evaluation and the
% function is not exited properly it could be possible, that the persistent
% variable is not deleted. In this case, call this function without any
% output arguments to delete the persistent variable and reset the display
% status.

% Persistent variables for output during evaluation
persistent opts;
persistent objfun;

% First function call from "bobyqa.m" => initialization
if (nargin == 3)
  
  % Objective function (handle)
  objfun = fun;%str2func(fun);
  
  % Set options
  cOpts = {'display', 'npt', 'rho_beg', 'rho_end', 'maxFunEval'};
  for i = 1:numel(cOpts)
    if isfield(options, cOpts{i})
      opts.(cOpts{i}) = options.(cOpts{i});
    else
      error(['Error setting the option structure for BOBYQA. The option ', cOpts{i}, 'is missing.']);
    end
  end
  
  % Initialize counter
  opts.nFunEval = 0;
  
% Final function call from "bobyqa.m" => reset svOpts
elseif (nargout == 0)
  clear svOpts;

% Function call from "bobyqa_alg.mex" => evaluate objective function and display status info
else
  fval = objfun(x);
  status_display()
end

%% Subfunctions
function status_display()
% STATUS_DISPLAY Displays the status in the Matlab command window depending
% on the 'display' option set in the options.

  switch opts.display
    case 'none'
    case 'iter'
      if mod(opts.nFunEval, 30) == 0
        fprintf(['\n', 'FunEval        ObjFunVal     Norm of step    Rel norm step', '\n']);
      end
      opts.nFunEval = opts.nFunEval + 1;
%       if isfield(opts, 'vX_last')
%         sNormStep = sprintf('%12.8e', norm(x - opts.vX_last));
%         sRelStep  = sprintf('%12.8e', norm(x - opts.vX_last)/opts.rho_end);
%       else
        sNormStep = '            ';
        sRelStep  = '            ';
%       end
      fprintf([sprintf('%7.1u', opts.nFunEval), '   ', ...
               sprintf('%12.8e', fval), '   ', ...
               sNormStep, '   ', ...
               sRelStep, ...
               '\n']);
      opts.vX_last = x;
    otherwise
      fprintf('Unknown display setting. Display set to ''none''.');
      opts.display = 'none';
  end
end

end