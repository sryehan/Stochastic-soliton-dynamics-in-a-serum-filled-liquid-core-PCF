function u = snlse_step_heun(u, k, dt, dB, kappa, c, varrho, gain)
%SNLSE_STEP_HEUN  One time step of  i L_t + (1/2) L_xx - kappa L - c |L|^2 L + varrho L o dB/dt = 0  (Stratonovich).
%   Strang splitting: half step of the dispersive term (exact in Fourier space), full step of the nonlinear,
%   detuning and noise terms with the stochastic Heun (predictor-corrector) scheme, half step of dispersion.
%   u      field, one row per realization (fft along dimension 2)
%   k      wave numbers (row vector), dt time step, dB Brownian increment (scalar or column vector, variance dt)
%   gain   optional linear gain; gain = varrho^2/2 reproduces the formal equation read in the Ito sense.
if nargin < 8, gain = 0; end
half = exp(-0.5i*k.^2*dt/2);
u = ifft(half.*fft(u, [], 2), [], 2);
f = @(v) (-1i*(kappa + c*abs(v).^2) + gain).*v;
g = @(v) 1i*varrho*v;
up = u + f(u)*dt + g(u).*dB;
u  = u + 0.5*(f(u) + f(up))*dt + 0.5*(g(u) + g(up)).*dB;
u = ifft(half.*fft(u, [], 2), [], 2);
end
