function dgmf = gmf_gradients(u, v, sigma, L, T, K, varargin)
% gmf_derivative defined the GMF gradients:
%
%	Fernando Cervantes Sanchez, June 2019
%	iie.fercer@gmail.com

    addpath('./GMF_functions/')
    delta = 1e-0;
    
    if numel(varargin) > 0
        delta = varargin{1};
    end
    
    dgmf.dsigma = cell(K, 1);
    dgmf.dL = cell(K, 1);
    dgmf.dT = cell(K, 1);
        
    for k = 0:(K-1)
        theta = k/K*pi;
        
        f = gmf_u_sigma(u, v, theta, sigma);
        df = dgmf_u_sigma(u, v, theta, sigma, f);
        
        h_v_L = step_fun_v_L(u, v, theta, L, delta);
        dh_v_L = dstep_fun_v_L(u, v, theta, L, delta);
                
        g_u_T = step_fun_u_T(u, v, theta, T, delta);
        dg_u_T = dstep_fun_u_T(u, v, theta, T, delta);
        
        active_volume = v_uv(u, v, theta, L, T, delta, h_v_L, g_u_T);
        dactive_volume = dv_uv(u, v, theta, L, T, delta, h_v_L, g_u_T, dh_v_L, dg_u_T);
        
        profile_volume = v_f_vu(u, v, theta, sigma, L, T, delta, f, h_v_L, g_u_T);
        dprofile_volume = dv_f_vu(u, v, theta, sigma, L, T, delta, f, h_v_L, g_u_T, df, dh_v_L, dg_u_T);
        
        dgmf.dsigma{k+1} = (((1-f)*dprofile_volume.dsigma - df.dsigma*(active_volume-profile_volume))/(active_volume-profile_volume)^2).*g_u_T.*h_v_L;
        
        dgmf.dL{k+1} = (-(1-f)*(dactive_volume.dL - dprofile_volume.dL)/(active_volume - profile_volume)^2 + dactive_volume.dL/active_volume^2).*g_u_T.*h_v_L ...
            + ((1-f)/(active_volume-profile_volume)-1/active_volume).*g_u_T.*dh_v_L.dL;
        
        dgmf.dT{k+1} = (-(1-f)*(dactive_volume.dT - dprofile_volume.dT)/(active_volume - profile_volume)^2 + dactive_volume.dT/active_volume^2).*g_u_T.*h_v_L ...
            + ((1-f)/(active_volume-profile_volume)-1/active_volume).*h_v_L.*dg_u_T.dT;        
    end
    rmpath('./GMF_functions/')
end