function [j] = sbesselj(nu,z)
% Spherical 

j = sqrt(pi/(2*z))*besselj(nu+0.5, z);

end
