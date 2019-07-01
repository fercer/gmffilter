function f = gmf_u_sigma(u, v, theta, sigma)
%gmf_u_sigma defines a gaussian profile along the x axis.
%
%   f = gmf_u_sigma(u, sigma) defines a Gaussian curve to model the profile
%   of a vessel in a medical image. The spread of the curve is defined by
%   sigma.
%
%   f = gmf_u_sigma(u, v, theta, sigma) defines a Gaussian curve to model 
%   the profile of a vessel in a medical image. The spread of the curve is 
%   defined by sigma. The curve is rotated theta radians.
    f = exp(-bsxfun(@rdivide, (u*cos(theta) - v*sin(theta)).^2,2*sigma.^2));
end