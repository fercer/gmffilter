function profile_volume = v_f_vu(u, v, theta, sigma, L, T, varargin)
%v_f_vu computes the volume under the profile of the Gaussian curve
%   delimited in the x and y axis.
%
%   profile_volume = v_f_vu(u, v, theta, sigma, L, T) computes the volume 
%   under the active profile of the Gaussian template.
%
%   profile_volume = v_f_vu(u, v, theta, sigma, L, T, delta, gmf_profile, g_u_T, h_v_L)
%   uses the precomputed profile of the Gaussian curve to save time. Also
%   the precomputed step functions are passes as arguments.
%
%   profile_volume = v_f_vu(u, v, theta, sigma, L, T, delta) computes the
%   the step functions using the delta value passed as argument.

    delta = 1e-0;
   
    if numel(varargin) == 1
        delta = varargin{1};
    end
    
    if numel(varargin) > 1
        f = varargin{2};
        h_v_L = varargin{3};
      	g_u_T = varargin{4};
    else
        f = gmf_u_sigma(u, v, theta, sigma);
        h_v_L = step_fun_v_L(u, v, theta, L, delta);
        g_u_T = step_fun_u_T(u, v, theta, T, delta);
    end
    
    U = u(1,:);
    V = v(:,1);
    
    profile_volume_fun = f.*g_u_T.*h_v_L;
    
    % Integral along u axis
    profile_volume_int_u = 0.5*sum(bsxfun(@times,profile_volume_fun(:,2:end,:)+profile_volume_fun(:,1:(end-1),:), diff(U)),2);
    
    profile_volume = 0.5*sum(bsxfun(@times,profile_volume_int_u(2:end,:)+profile_volume_int_u(1:(end-1),:), diff(V)));
end