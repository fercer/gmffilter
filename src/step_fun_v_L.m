function h = step_fun_v_L(u, v, theta, L, varargin)
%step_fun_v_L defines a smooth differentiable step function traslated L 
%   units up on the y axis.
%
%   h = step_fun_v_L(u, v, theta, L) defines a step function to delimit a 
%   function on the y axis (rotated theta radians). A delta of 1e-0 is 
%   assumed if nothing is passed.
%
%   h = step_fun_v_L(u, v, theta, L, delta) defines a step function to 
%   delimit a function on the y axis (rotated theta radians). The delta 
%   passed by the user defines how smooth or stepest is the function.
%
%  This version of step_fun_v_L used the logistic sigmoid function as the 
%  step function. It can be changed for any other differentiable step 
%  function.

    delta = 1e-0;
       
    if numel(varargin) > 0
        delta = varargin{1};
    end

    h = 1./(1 + exp(-bsxfun(@plus, v*cos(theta)+u*sin(theta), L/2)/delta)) - ...
        1./(1 + exp(-bsxfun(@minus, v*cos(theta)+u*sin(theta), L/2)/delta));
end