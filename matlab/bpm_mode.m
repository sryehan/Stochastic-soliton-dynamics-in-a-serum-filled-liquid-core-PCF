function [psi, neff, corefrac, Aeff] = bpm_mode(eps, x, lam0, n0, rmask, w0, Rcore)
%BPM_MODE  Fundamental core mode of a cross-section by beam propagation in imaginary distance (z -> -i*tau).
%   In imaginary distance every mode decays as exp(-(k0 n0 - beta) tau), so a Gaussian launched in the core
%   relaxes to the mode with the largest propagation constant beta.  The silica beyond the radius rmask is
%   replaced by air for this search; otherwise the solid jacket (index 1.45), which has the highest index of
%   the whole cross-section, would win.  The core mode does not reach that far, so it is not affected.
%
%   psi       mode field, real, normalised to sum(psi(:).^2) = 1
%   neff      effective index, from the expectation value of the operator (scalar Helmholtz value)
%   corefrac  fraction of the mode power inside r <= Rcore
%   Aeff      effective mode area (um^2), (sum I)^2 / sum I^2 * dx^2
N = numel(x); dx = x(2) - x(1); k0 = 2*pi/lam0;
[X, Y] = meshgrid(x, x); R = hypot(X, Y);
eps(R > rmask) = 1.0;
kx = 2*pi/(N*dx)*[0:N/2-1, -N/2:-1]; [KX, KY] = meshgrid(kx, kx);
T = (KX.^2 + KY.^2)/(2*k0*n0);                                   % "kinetic" part (1/um)
U = -k0*(eps - n0^2)/(2*n0);                                     % "potential" part (1/um)
psi = exp(-R.^2/w0^2);
dtau = [0.2, 0.05, 0.0125]; nst = [400, 600, 800];               % coarse to fine: total imaginary distance 120 um
for s = 1:numel(dtau)
    eT = exp(-T*dtau(s)); eU = exp(-U*dtau(s)/2);
    for it = 1:nst(s)
        psi = real(eU.*ifft2(eT.*fft2(eU.*psi)));
        psi = psi/sqrt(sum(psi(:).^2));
    end
end
F = fft2(psi);
E = real(sum(conj(F(:)).*T(:).*F(:))/sum(abs(F(:)).^2)) + sum(U(:).*psi(:).^2);   % = k0 n0 - beta (paraxial)
neff = sqrt(n0^2 - 2*n0*E/k0);
I = psi.^2;
corefrac = sum(I(R <= Rcore))/sum(I(:));
Aeff = sum(I(:))^2/sum(I(:).^2)*dx^2;
end
