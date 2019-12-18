function dactive_volume = dv_uv(u, v, theta, L, T, varargin)
%dv_uv computes the gradient of the volume under the difference of two
%   step functions, on two different axes x and y.
%
%   vol_profile = dv_uv(u, v, theta, L, T)  computes the gradient of the 
%   volume between the step functions defined by T on the x axis, and L in
%   the y axis. The axes are rotated theta radians.
%
%   vol_profile = dv_uv(u, v, ..., delta)  computes the gradient of the 
%   volume between the step functions defined by T on the x axis, and L in
%   the y axis using the delta defined by the user.
%
%   vol_profile = dv_uv(u, v, ..., h_v_L, g_u_T, dh_v_L, dg_u_T) computes
%   the gradient of the volume using the precomputed functions and
%   gradients of the step functions.
    
    delta = 1e-0;
   
    if numel(varargin) == 1
        delta = varargin{1};
    end
    
    if numel(varargin) > 1
        h_v_L = varargin{2};
        g_u_T = varargin{3};
        
        dh_v_L = varargin{4};
        dg_u_T = varargin{5};
    else    
        % Compute the setep function for the rotated y axis:
        h_v_L = step_fun_v_L(u, v, theta, L, delta);
        dh_v_L = dstep_fun_v_L(u, v, theta, L, delta);
        
        % Compute the setep function for the rotated x axis:
        g_u_T = step_fun_u_T(u, v, theta, T, delta);
        dg_u_T = dstep_fun_u_T(u, v, theta, T, delta);
    end
    
    dvolume_grdL = g_u_T.*dh_v_L.dL;
    dvolume_grdT = h_v_L.*dg_u_T.dT;
    
    U = u(1,:);
    V = v(:,1);
    
    % Integral along u axis
    dactive_volume_grdL_int_u = 0.5*sum(bsxfun(@times,dvolume_grdL(:,2:end,:)+dvolume_grdL(:,1:(end-1),:), diff(U)),2);
    
    % Integral along v axis
    dactive_volume.dL = 0.5*sum(bsxfun(@times,dactive_volume_grdL_int_u(2:end,:)+dactive_volume_grdL_int_u(1:(end-1),:), diff(V)));
    
    % Integral along u axis    
    dactive_volume_grdT_int_u = 0.5*sum(bsxfun(@times,dvolume_grdT(:,2:end,:)+dvolume_grdT(:,1:(end-1),:), diff(U)),2);
    
    % Integral along v axis
    dactive_volume.dT = 0.5*sum(bsxfun(@times,dactive_volume_grdT_int_u(2:end,:)+dactive_volume_grdT_int_u(1:(end-1),:), diff(V))); 
end