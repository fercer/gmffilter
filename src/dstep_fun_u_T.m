function dg = dstep_fun_u_T(u, v, theta, T, varargin)
%dstep_fun_u_T defines the gradient of the step function step_fun_u_T.
%
%   dg = dstep_fun_u_T(u, v, theta, T) defines the gradient of a step
%   function that delimits a function on the x axis (rotated theta radians).
%   A delta of 1e-0 is assumed if nothing is passed.
%
%   dg = dstep_fun_u_T(u, v, theta, T, delta)  defines the gradient of
%   a step function that delimits a function on the x axis (rotated theta 
%   radians). The delta passed by the user defines how smooth or stepest 
%   is the function.
%
%	This version of dstep_fun_u_T used the logistic sigmoid function as the
%   step function. It can be changed for any other differentiable step
%   function.

    delta = 1e-0;
    
    if numel(varargin) > 0
        delta = varargin{1};
    end
    
    g_inf = 1./(1 + exp(-(u*cos(theta) - v*sin(theta) + T/2)/delta));
    g_sup = 1./(1 + exp(-(u*cos(theta) - v*sin(theta) - T/2)/delta));
    
    dg_inf = g_inf.*(1-g_inf)/(2*delta);
    dg_sup = -g_sup.*(1-g_sup)/(2*delta);
    
    dg.dT = dg_inf - dg_sup;
end