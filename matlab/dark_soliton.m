function L = dark_soliton(x, t, Bt, K, m, kappa, c, varrho, x0, l)
%DARK_SOLITON  Exact dark soliton of the stochastic NLSE (Eq. 27 of the paper), c > 0, K = k*w.
if nargin < 9, x0 = 0; end
if nargin < 10, l = 0; end
n = -0.5*(2*K^2 + m^2 + 2*kappa);
L = K/sqrt(c)*tanh(K*(x - x0 - m*t)).*exp(1i*(m*x + n*t + l + varrho*Bt));
end
