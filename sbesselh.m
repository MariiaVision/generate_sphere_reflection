function [j] = sbesselh(nu,z)
% Spherical 

j = sqrt(pi/(2*z))*besselh(nu+0.5, z);

end
