function [j] = dbesselh(nu,z)
%DBESSELH First derivative of the Henkel function of the first kind.
% spherical case

j = 0.5*(sbesselh(nu-1,z) - sbesselh(nu+1,z));

end
