function [neff, alphaP] = capillary_leaky_mode(lam0, a, n_core, n_clad)
%CAPILLARY_LEAKY_MODE  Fundamental leaky mode of a round core of radius a whose index is below that of the cladding.
%   Scalar wave equation.  The field is J0(u r/a) in the core and the outgoing Hankel function H0(w r/a) outside,
%   u = a sqrt(k0^2 n_core^2 - beta^2),  w = a sqrt(k0^2 n_clad^2 - beta^2).  Continuity of the field and of its
%   radial derivative gives   u J1(u) H0(w) = w H1(w) J0(u),   solved here for the complex propagation constant
%   beta by the secant method, starting from the estimate of Marcatili and Schmeltzer (Bell Syst. Tech. J. 43,
%   1783, 1964).
%   neff = Re(beta)/k0;  alphaP = 2 Im(beta) = decay rate of the power, per unit of length (unit of lam0 and a).
k0 = 2*pi/lam0; nu = n_clad/n_core; lamc = lam0/n_core; u0 = 2.405;
U = @(b) a*sqrt(k0^2*n_core^2 - b.^2);
W = @(b) a*sqrt(k0^2*n_clad^2 - b.^2);
F = @(b) U(b).*besselj(1, U(b)).*besselh(0, 1, W(b)) - W(b).*besselh(1, 1, W(b)).*besselj(0, U(b));
b0 = k0*n_core*(1 - 0.5*(u0*lamc/(2*pi*a))^2) + 1i*(u0/(2*pi))^2*lamc^2/a^3/sqrt(nu^2 - 1);
b1 = b0*(1 + 1e-5); f0 = F(b0); f1 = F(b1);
for it = 1:60
    b2 = b1 - f1*(b1 - b0)/(f1 - f0);
    b0 = b1; f0 = f1; b1 = b2; f1 = F(b1);
    if abs(b1 - b0) < 1e-13*abs(b1), break; end
end
neff = real(b1)/k0; alphaP = 2*imag(b1);
end
