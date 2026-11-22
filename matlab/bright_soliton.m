function L = bright_soliton(x, t, Bt, K, m, kappa, c, varrho, x0, l)
%BRIGHT_SOLITON  Exact bright soliton of the stochastic NLSE (Eq. 26 of the paper), c < 0, K = k*w.
%   Bt is the value of the Brownian path B(t) (0 for the deterministic solution).
if nargin < 9, x0 = 0; end
if nargin < 10, l = 0; end
n = 0.5*(K^2 - m^2 - 2*kappa);
L = K/sqrt(-c)*sech(K*(x - x0 - m*t)).*exp(1i*(m*x + n*t + l + varrho*Bt));
end
