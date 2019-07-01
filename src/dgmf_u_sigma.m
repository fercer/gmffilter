function df = dgmf_u_sigma(u, v, theta, sigma, varargin)
%dgmf_u_sigma defines the gradients of the gaussian profile along the x axis.
%
%   df = dgmf_u_sigma(u, v, theta, sigma) defines the gradients of the Gaussian
%   curve rotated theta degrees.
%
%   df = dgmf_u_sigma(u, v, theta, sigma, varargin) defines the gradients
%   of the Gaussian curve rotated theta degrees. A precomputed Gaussian
%   profile is used as passed by the user.

    if numel(varargin) > 0
        f = varargin{1};
    else
        f = gmf_u_sigma(u, v, theta, sigma);
    end
    
    df.dsigma = (u*cos(theta) - v*sin(theta)).^2.*f/sigma^3;
end