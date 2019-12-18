function g = step_fun_u_T(u, v, theta, T, varargin)
%step_fun_u_T defines a smooth differentiable step function traslated T units
%   to the right on the x axis.
%
%   g = step_fun_u_T(u, theta, T) defines a step function to delimit a 
%   function on the x axis (rotated theta radians). A delta of 1e-0 is 
%   assumed if nothing is passed.
%
%   g = step_fun_u_T(u, theta, T, delta) defines a step function to 
%   delimit a function on the x axis (rotated theta radians). The delta 
%   passed by the user defines how smooth or stepest is the function.
%
%	This version of step_fun_u_T used the logistic sigmoid function as the
%   step function. It can be changed for any other differentiable step 
%   function.

    delta = 1e-0;
    
    if numel(varargin) > 0
        delta = varargin{1};
    end

    g = 1./(1 + exp(-bsxfun(@plus, u*cos(theta) - v*sin(theta), T/2)/delta)) - ...
        1./(1 + exp(-bsxfun(@minus, u*cos(theta) - v*sin(theta), T/2)/delta));
end