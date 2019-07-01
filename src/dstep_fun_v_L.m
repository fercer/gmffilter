function dh = dstep_fun_v_L(u, v, theta, L, varargin)
%dstep_fun_v_L defines the gradient of the step function step_fun_v_TL
%
%   dh = dstep_fun_v_L(u, T) defines the gradient of a step function 
%   that delimits a function on the y axis. A delta of 1e-0 is assumed if 
%   nothing is passed.
%
%   dh = dstep_fun_v_L(u, v, theta, L) defines the gradient of a step
%   function that delimits a function on the y axis (rotated theta radians).
%   A delta of 1e-12 is assumed if nothing is passed.
%
%   dh = dstep_fun_v_L(u, v, theta, L, delta)  defines the gradient of
%   a step function that delimits a function on the y axis (rotated theta 
%   radians). The delta passed by the user defines how smooth or stepest 
%   is the function.
%
%	This version of dstep_fun_v_L used the logistic sigmoid function as the
%   step function. It can be changed for any other differentiable step
%   function.

    delta = 1e-0;    
   
    if numel(varargin) > 0
        delta = varargin{1};
    end
    
    h_inf = 1./(1 + exp(-(v*cos(theta) + u*sin(theta) + L/2)/delta));
    h_sup = 1./(1 + exp(-(v*cos(theta) + u*sin(theta) - L/2)/delta));
    
    dh_inf = h_inf.*(1-h_inf)/(2*delta);
    dh_sup = -h_sup.*(1-h_sup)/(2*delta);
    
    dh.dL = dh_inf - dh_sup;
end