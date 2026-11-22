function out = bpm_run(eps, x, lam0, n0, dz, Z, psi0, Rcore, rabs, zsnap)
%BPM_RUN  Scalar paraxial beam propagation (split-step Fourier) through a z-invariant fiber cross-section.
%   Solves   d(psi)/dz = i/(2 k0 n0) * (d2/dx2 + d2/dy2) psi + i k0 (n^2 - n0^2)/(2 n0) * psi ,
%   the paraxial wave equation for the field E = psi * exp(i k0 n0 z).  It is a two-dimensional Schroedinger
%   equation in which the distance z plays the role of time (z = c0 t / n_g for a pulse of group index n_g).
%
%   eps   n^2(x,y) on the square grid x (um), e.g. from pcf_index_map
%   lam0  vacuum wavelength (um);  n0  reference index;  dz  step (um);  Z  length (um)
%   psi0  launched field: either an N x N array or a scalar = 1/e field radius (um) of a centred Gaussian
%   Rcore radius (um) used for the "power in the core";  rabs  radius at which the absorbing boundary starts
%   zsnap distances (um) at which the intensity |psi|^2 is stored (optional)
%
%   out.z      distances (um), one record per micrometre
%   out.Pcore  power inside r <= Rcore, as a fraction of the launched power
%   out.Ptot   power left in the whole window (the absorber removes what has radiated away)
%   out.c      overlap <psi0|psi(z)>/<psi0|psi0>; its phase advances at the rate k0 (n_eff - n0)
%   out.snap   cell array of intensity maps at out.snap_z, normalised to the launched peak intensity
%
%   The step must be small because the index jumps from 1.45 (silica) to 1.00 (air): with dz = 0.2 um the
%   splitting error scatters about 9 % of the power per 100 um out of design C, with dz = 0.025 um 0.1-0.2 %.
%   The absorber rises quadratically from zero at rabs to 1 per um at the edge of the window.  Keep it about
%   5 um thick: a thin, strong absorber reflects part of the radiated light back into the fiber.
if nargin < 10, zsnap = []; end
N = numel(x); dx = x(2) - x(1); k0 = 2*pi/lam0;
[X, Y] = meshgrid(x, x); R = hypot(X, Y);
kx = 2*pi/(N*dx)*[0:N/2-1, -N/2:-1]; [KX, KY] = meshgrid(kx, kx);
D  = exp(-1i*(KX.^2 + KY.^2)/(2*k0*n0)*dz);                      % diffraction, full step
Dh = exp(-1i*(KX.^2 + KY.^2)/(2*k0*n0)*dz/2);                    % diffraction, half step
rmax = max(x); sig = 1.0*max(0, (R - rabs)/(rmax - rabs)).^2;    % absorber strength (1/um), zero inside rabs
V = exp(1i*k0*(eps - n0^2)/(2*n0)*dz - sig*dz);                  % index profile and absorber, full step
if isscalar(psi0), psi0 = exp(-R.^2/psi0^2); end
P0 = sum(abs(psi0(:)).^2); I0 = max(abs(psi0(:)).^2);
core = R <= Rcore;
nst = round(Z/dz); every = max(1, round(1/dz)); nrec = floor(nst/every);
out.z = zeros(1, nrec); out.Pcore = zeros(1, nrec); out.Ptot = zeros(1, nrec); out.c = zeros(1, nrec);
isnap = round(zsnap/dz); out.snap = cell(1, numel(zsnap)); out.snap_z = isnap*dz;
psi = ifft2(Dh.*fft2(psi0));                                     % first half step of diffraction
ir = 0;
for it = 1:nst
    psi = psi.*V;
    rec = mod(it, every) == 0; q = find(isnap == it);
    if rec || ~isempty(q)
        f = ifft2(Dh.*fft2(psi));                                % field at z = it*dz
        I = abs(f).^2;
        if rec
            ir = ir + 1; out.z(ir) = it*dz;
            out.Pcore(ir) = sum(I(core))/P0; out.Ptot(ir) = sum(I(:))/P0;
            out.c(ir) = sum(conj(psi0(:)).*f(:))/P0;
        end
        for j = q, out.snap{j} = I/I0; end
    end
    psi = ifft2(D.*fft2(psi));
end
end
