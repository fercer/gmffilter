function active_volume = v_uv(u, v, theta, L, T, varargin)
%v_uv computes the volume under the difference of two step functions, on
%   two different axes x and y.    
%
%   vol_profile = v_uv(u, v, theta, L, T) computes the volume between 
%   the step functions defined by T on the x axis, and L in the y axis. The
%   axes are rotated theta radians.
%
%   vol_profile = v_uv(u, v, ..., delta) computes the volume between the 
%   step functions defined by T on the x axis, and L in the y axis using 
%   the delta defined by the user.
%
%   vol_profile = v_uv(u, v, ..., delta) computes the volume using the
%   precomputed step functions to save time.
    
    delta = 1e-0;
   
    if numel(varargin) == 1
        delta = varargin{1};
    end
    
    if numel(varargin) > 1
        h_v_L = varargin{2};
        g_u_T = varargin{3};
    else
        h_v_L = step_fun_v_L(u, v, theta, L, delta);
        g_u_T = step_fun_u_T(u, v, theta, T, delta);
    end
    
    U = u(1,:);
    V = v(:,1);
    
    
    active_volume_fun = g_u_T.*h_v_L;
    
    % Integral along u axis
    actuve_volume_int_u = 0.5*sum(bsxfun(@times, active_volume_fun(:,2:end,:)+active_volume_fun(:,1:(end-1),:), diff(U)), 2);
    
    %Integral along v axis
    active_volume = 0.5*sum(bsxfun(@times, actuve_volume_int_u(2:end,:)+actuve_volume_int_u(1:(end-1),:), diff(V)),1);
end