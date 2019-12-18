function dprofile_volume = dv_f_vu(u, v, theta, sigma, L, T, varargin)
%dv_f_vu computes the gradient of the volume under the profile of the 
%   Gaussian curve delimited in the x and y axis.
%
%   dprofile_volume = dv_f_vu(u, v, theta, sigma, L, T) computes the 
%   gradient volume under the active profile of the Gaussian template.
%
%   dprofile_volume = dv_f_vu(u, v, theta, sigma, L, T, delta, gmf_profile, g_u_T, h_v_L, dgmf_profile, dg_u_T, dh_v_L)
%   uses the precomputed profile of the Gaussian curve to save time. Also
%   the precomputed step functions are passes as arguments and all their
%   gradients can be passes as arguments.
%
%   dprofile_volume = dv_f_vu(u, v, theta, sigma, L, T, delta) computes
%   the gradient of the area under the profile delimited by the step
%   functions computed using the delta value passed as argument.

    delta = 1e-0;
   
    if numel(varargin) == 1
        delta = varargin{1};
    end
    
    if numel(varargin) > 1
        f = varargin{2};
        h_v_L = varargin{3};
      	g_u_T = varargin{4};
        
        df = varargin{5};
        dh_v_L = varargin{6};
      	dg_u_T = varargin{7};
    else
        f = gmf_u_sigma(u, v, theta, sigma);
        df = dgmf_u_sigma(u, v, theta, sigma, f);
            
        % Compute the setep function for the rotated y axis:
        h_v_L = step_fun_v_L(u, v, theta, L, delta);
        dh_v_L = dstep_fun_v_L(u, v, theta, L, delta);
        
        % Compute the setep function for the rotated x axis:
        g_u_T = step_fun_u_T(u, v, theta, T, delta);
        dg_u_T = dstep_fun_u_T(u, v, theta, T, delta);
    end
            
    % Compute the setep function for the rotated x axis:
    U = u(1,:);
    V = v(:,1);
    
    dprofile_volume_grdsigma = h_v_L.*g_u_T.*df.dsigma;
    dprofile_volume_grdL = f.*g_u_T.*dh_v_L.dL;
    dprofile_volume_grdT = f.*h_v_L.*dg_u_T.dT;
    
    % Integral along u axis
    dprofile_volume_grdsigma_int_u = 0.5*sum(bsxfun(@times,dprofile_volume_grdsigma(:,2:end,:)+dprofile_volume_grdsigma(:,1:(end-1),:), diff(U)),2);
    
    % Integral along v axis
    dprofile_volume.dsigma = 0.5*sum(bsxfun(@times,dprofile_volume_grdsigma_int_u(2:end,:)+dprofile_volume_grdsigma_int_u(1:(end-1),:), diff(V)));
    
    % Integral along u axis
    dprofile_volume_grdL_int_u = 0.5*sum(bsxfun(@times,dprofile_volume_grdL(:,2:end,:)+dprofile_volume_grdL(:,1:(end-1),:), diff(U)),2);
    
    % Integral along v axis
    dprofile_volume.dL = 0.5*sum(bsxfun(@times,dprofile_volume_grdL_int_u(2:end,:)+dprofile_volume_grdL_int_u(1:(end-1),:), diff(V)));
    
    % Integral along u axis
    dprofile_volume_grdT_int_u = 0.5*sum(bsxfun(@times,dprofile_volume_grdT(:,2:end,:)+dprofile_volume_grdT(:,1:(end-1),:), diff(U)),2);
    
    % Integral along v axis
    dprofile_volume.dT = 0.5*sum(bsxfun(@times,dprofile_volume_grdT_int_u(2:end,:)+dprofile_volume_grdT_int_u(1:(end-1),:), diff(V)));    
end