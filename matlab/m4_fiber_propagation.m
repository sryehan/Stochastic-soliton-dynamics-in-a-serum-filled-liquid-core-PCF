% M4_FIBER_PROPAGATION  Fundamental soliton in the serum-filled liquid-core PCF with absorption and third-order
% dispersion, compared with linear propagation of the same pulse (Section 6.4, Fig. 8).  Split-step Fourier method.
% Modal parameters at 1064 nm from the mode solver (out/numbers.json).
clear; close all;
beta2 = -26.6081e-27;      % s^2/m   (-26.6081 ps^2/km)
beta3 = 0.13535e-39;      % s^3/m   (0.13535 ps^3/km)
gamma = 0.020548;          % 1/(W m)
alpha = 14.1484;          % 1/m  (absorption of water times the index sensitivity)
for T0 = [30e-15 20e-15]
    LD = T0^2/abs(beta2); P0 = abs(beta2)/(gamma*T0^2); Gam = alpha*LD/2; d3 = beta3/(6*abs(beta2)*T0);
    fprintf('T0 = %2.0f fs: FWHM = %.1f fs, L_D = %.3f cm, P0 = %.1f W, energy = %.1f pJ, Gamma = %.3f, delta3 = %.4f\n', ...
            T0*1e15, 1.7627*T0*1e15, LD*100, P0, 2*P0*T0*1e12, Gam, d3);
    N = 4096; Lx = 200; dxi = 5e-4; x = (-N/2:N/2-1)*(Lx/N); k = 2*pi/Lx*[0:N/2-1, -N/2:-1];
    half = exp((-Gam - 0.5i*k.^2 - 1i*d3*k.^3)*dxi/2);          % loss, dispersion and TOD, half step
    targets = [2 5 8]/100/LD; nsteps = round(max(targets)/dxi);
    us = sech(x); ul = us; f0 = pulse_fwhm(x, abs(us).^2); E0 = sum(abs(us).^2);
    for it = 1:nsteps
        us = ifft(half.*fft(us)); us = us.*exp(1i*abs(us).^2*dxi); us = ifft(half.*fft(us));   % soliton (N = 1)
        ul = ifft(half.*fft(ul)); ul = ifft(half.*fft(ul));                                    % linear reference
        j = find(abs(targets - it*dxi) < dxi/2);
        if ~isempty(j)
            fprintf('   z = %d cm: FWHM/FWHM0 soliton %.3f, linear %.3f; energy %.3f, peak power %.3f\n', ...
                    round(targets(j)*LD*100), pulse_fwhm(x, abs(us).^2)/f0, pulse_fwhm(x, abs(ul).^2)/f0, sum(abs(us).^2)/E0, max(abs(us).^2));
        end
    end
end
