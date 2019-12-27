function gmf = gmf_functions(u, v, sigma, L, T, K, varargin)
% gmf_functions defined the GMF functions required to filter an image:
%
%	Fernando Cervantes Sanchez, June 2019
%	iie.fercer@gmail.com
    
    addpath('./GMF_functions/')
    delta = 1e-0;
    
    if numel(varargin) > 0
        delta = varargin{1};
    end
        
    [rows, cols] = size(u);
    
    gmf = cell(K,1);
    for k = 0:(K-1)
        theta = k/K*pi;
        
        f = gmf_u_sigma(u, v, theta, sigma);        
        h_v_L = step_fun_v_L(u, v, theta, L, delta);
        g_u_T = step_fun_u_T(u, v, theta, T, delta);
        
        active_volume = v_uv(u, v, theta, L, T, delta, h_v_L, g_u_T);
        profile_volume = v_f_vu(u, v, theta, sigma, L, T, delta, f, h_v_L, g_u_T);
        
        active_volume = repmat(active_volume, [rows, 1, cols]);
        active_volume = permute(active_volume,[1,3,2]);
        
        profile_volume = repmat(profile_volume, [rows, 1, cols]);
        profile_volume = permute(profile_volume,[1,3,2]);
        
        gmf{k+1} = ((1-f)./(active_volume-profile_volume) - 1./active_volume).*g_u_T.*h_v_L;
    end
    
    rmpath('./GMF_functions/')
end