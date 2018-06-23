function [j] = ddbessely(nu,z)
%DDBESSELJ Second derivative of the Bessel function of the second kind.

% spherical case

j = 0.25*(sbessely(nu-2,z) - 2*sbessely(nu,z)+sbessely(nu+2,z));

end
