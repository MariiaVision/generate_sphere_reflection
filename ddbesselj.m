function [j] = ddbesselj(nu,z)
%DDBESSELJ Second derivative of the Bessel function of the second kind.

% spherical case

j = 0.25*(sbesselj(nu-2,z) - 2*sbesselj(nu,z)+sbesselj(nu+2,z));

end
