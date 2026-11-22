% M6_BEAM_PROPAGATION  Light guidance over distance in three liquid-core fibers, by the beam propagation method.
%
% The same Gaussian beam (1064 nm) is launched into the liquid core of
%   A  a liquid-filled capillary: liquid core in solid silica,
%   B  a conventional photonic crystal fiber cladding (pitch 4.46 um, d/Lambda = 0.5) with a liquid core,
%   C  the fiber of this work (pitch 1.66 um, d/Lambda = 0.952, 80 nm core wall, 120 air holes),
% and the paraxial wave equation
%   d(psi)/dz = i/(2 k0 n0) (d2/dx2 + d2/dy2) psi + i k0 (n(x,y)^2 - n0^2)/(2 n0) psi
% is integrated along the fiber with the split-step Fourier method (bpm_run.m).  The distance z plays the role
% of time; a pulse needs t = z n_g/c0 to travel the distance z.
%
% The script then finds the fundamental mode of design C (bpm_mode.m), compares it with the full-vector mode
% solver of the paper, and propagates that mode through design C filled with 80 and 300 mg/dL glucose to show
% how the optical phase difference builds up along the fiber.
%
% The model is scalar and paraxial.  It shows whether a structure guides, and it is an independent check of the
% design; the numbers of the paper come from the full-vector mode solver, not from this script.
% Run time: 9 min in GNU Octave 8.4 on a two-core machine (fast = true: about half).  Written for MATLAB and
% Octave with base functions only; it has been executed in Octave (output in octave_output_m6.txt).
% Needs pcf_index_map.m, bpm_run.m, bpm_mode.m and capillary_leaky_mode.m in the same folder.
clear; close all;
fast = false;                                 % true: step 0.05 um instead of 0.025 um (about 1 % numerical loss per 100 um)

% ---------------- parameters ----------------
lam0 = 1.064; k0 = 2*pi/lam0;                 % wavelength (um)
n_si = 1.449631;                              % fused silica at 1064 nm (Malitson)
n_w  = 1.323917;                              % water at 1064 nm (Kedenburg et al.)
n_a  = n_w + 0.016;                           % serum-like liquid at 100 mg/dL glucose
dndC = 1.4e-6; Cref = 100;                    % index change per mg/dL of glucose; reference concentration
n0 = 1.332;                                   % reference index of the paraxial equation
N = 640; dx = 0.05; x = (-N/2:N/2-1)*dx;      % 32 um window, 50 nm pixels
dz = 0.025; if fast, dz = 0.05; end           % propagation step (um)
Z = 200;                                      % fiber length for the guidance test (um), at least 100
Zg = 100;                                     % fiber length for the glucose test (um)
w0 = 1.9;                                     % 1/e field radius of the launched Gaussian (um)
Rc = 2.45;                                    % radius of the liquid core (um)
rabs = 11.0;                                  % the absorbing boundary starts at this radius (um), 5 um thick
zsnap = [20 50 200];                          % distances of the intensity pictures (um)
c0 = 299.792458; n_g = 1.3633;                % speed of light (um/ps); group index of the mode of design C
% full-vector mode solver of the paper (design C, out/numbers.json), for comparison
vec.neff = 1.332018; vec.core = 0.9919; vec.Aeff = 11.78; vec.Sn = 0.99638; vec.dphi = 90.61;

names = {'A: liquid core in solid silica', 'B: conventional cladding, d/\Lambda = 0.5', 'C: this work, d/\Lambda = 0.952'};
des = 'ABC'; maps = cell(1, 3); res = cell(1, 3);

% ---------------- part 1: the same beam through the three fibers ----------------
fprintf('grid %d x %d, pixel %.0f nm, step %.0f nm, length %.0f um = %.2f ps of travel time\n', N, N, 1e3*dx, 1e3*dz, Z, Z*n_g/c0);
for q = 1:3
    maps{q} = pcf_index_map(des(q), x, n_si, n_a);
    tic; res{q} = bpm_run(maps{q}, x, lam0, n0, dz, Z, w0, Rc, rabs, zsnap);
    fprintf('design %s done in %.0f s\n', des(q), toc);
end
fprintf('\npower inside the liquid core / power left in the window (fractions of the launched power)\n');
fprintf('   z (um)        A                  B                  C\n');
zlist = [10 20 50 100 150 200]; zlist = zlist(zlist <= Z);
for zq = zlist
    j = find(res{1}.z >= zq - 1e-9, 1);
    fprintf('   %5.0f    %.4f / %.4f    %.4f / %.4f    %.4f / %.4f\n', zq, res{1}.Pcore(j), res{1}.Ptot(j), ...
            res{2}.Pcore(j), res{2}.Ptot(j), res{3}.Pcore(j), res{3}.Ptot(j));
end
% leakage rates: straight-line fits of ln(power) between 50 um and the end
dB = 10/log(10)*1e3;                          % from 1/um to dB/mm
m = res{1}.z >= 50; pA = polyfit(res{1}.z(m), log(res{1}.Pcore(m)), 1); aA = -pA(1);
m = res{2}.z >= 50; pB = polyfit(res{2}.z(m), log(res{2}.Ptot(m)), 1);  aB = -pB(1);
[nA, aEx] = capillary_leaky_mode(lam0, Rc, n_a, n_si);              % exact leaky mode of the capillary (scalar)
nu = n_si/n_a; lamc = lam0/n_a;
aMS_s = 2*(2.405/(2*pi))^2*lamc^2/Rc^3/sqrt(nu^2 - 1);            % Marcatili-Schmeltzer estimate, scalar
aMS_v = aMS_s*(nu^2 + 1)/2;                                       % same, hybrid EH11 mode of the vector problem
fprintf('\nA: core power decays at %.4f per um = %.0f dB/mm (1/e length %.0f um)\n', aA, aA*dB, 1/aA);
fprintf('   exact leaky mode of the capillary: %.4f per um; Marcatili-Schmeltzer estimate: %.4f (scalar), %.4f (EH11) per um\n', aEx, aMS_s, aMS_v);
fprintf('B: power in the window decays at %.4f per um = %.0f dB/mm on average; core power %.3f at 50 um, %.3f at %.0f um\n', ...
        aB, aB*dB, res{2}.Pcore(find(res{2}.z >= 50, 1)), res{2}.Pcore(end), Z);
m = res{3}.z >= Z/2;
fprintf('C: core power %.4f +- %.4f over the last %.0f um; power left in the window at %.0f um: %.4f\n', ...
        mean(res{3}.Pcore(m)), std(res{3}.Pcore(m)), Z/2, Z, res{3}.Ptot(end));

% ---------------- part 2: fundamental mode of design C ----------------
tic; [psiC, neffC, coreC, AeffC] = bpm_mode(maps{3}, x, lam0, n0, 9.4, w0, Rc);
[X, Y] = meshgrid(x, x); g = exp(-(X.^2 + Y.^2)/w0^2);
eta = sum(g(:).*psiC(:))^2/sum(g(:).^2);
fprintf('\nfundamental mode of design C (%.0f s)            this script (scalar)    paper (full-vector)\n', toc);
fprintf('   effective index                               %.6f                %.6f\n', neffC, vec.neff);
fprintf('   power fraction in the liquid                  %.4f                  %.4f\n', coreC, vec.core);
fprintf('   effective area (um^2)                         %.2f                   %.2f\n', AeffC, vec.Aeff);
fprintf('   difference of the effective index: %.1e;  overlap of the launched Gaussian with the mode: %.4f\n', neffC - vec.neff, eta);

% ---------------- part 3: phase difference between two glucose levels in design C ----------------
Cg = [80 300]; og = cell(1, 2);
for q = 1:2
    eg = pcf_index_map('C', x, n_si, n_a + dndC*(Cg(q) - Cref));
    tic; og{q} = bpm_run(eg, x, lam0, n0, dz, Zg, psiC, Rc, rabs);
    fprintf('glucose %d mg/dL done in %.0f s\n', Cg(q), toc);
end
zg = og{1}.z; dphi = unwrap(angle(og{2}.c.*conj(og{1}.c)));
pg = polyfit(zg, dphi, 1);
ph1 = unwrap(angle(og{1}.c)); p1 = polyfit(zg, ph1, 1);
Sn = pg(1)/(k0*dndC*(Cg(2) - Cg(1)));
fprintf('\neffective index from the phase advance at %d mg/dL: %.6f\n', Cg(1), n0 + p1(1)/k0);
fprintf('phase difference between %d and %d mg/dL: %.3e rad per um, i.e. %.4f rad after %.0f um\n', Cg(1), Cg(2), pg(1), pg(1)*Zg, Zg);
fprintf('   d n_eff/d n_a = %.4f (paper %.4f);  extrapolated to 5 cm: %.1f rad = %.1f fringes (paper %.1f rad)\n', ...
        Sn, vec.Sn, pg(1)*5e4, pg(1)*5e4/(2*pi), vec.dphi);
fprintf('   power left after %.0f um: %.4f (80 mg/dL), %.4f (300 mg/dL)\n', Zg, og{1}.Ptot(end), og{2}.Ptot(end));

save('bpm_results.mat', 'x', 'maps', 'res', 'psiC', 'neffC', 'coreC', 'AeffC', 'eta', 'aA', 'aB', 'aEx', 'aMS_s', 'aMS_v', ...
     'zg', 'dphi', 'pg', 'Sn', 'og', 'Cg', 'Z', 'Zg', 'dz', 'w0', '-v7');

% ---------------- figures ----------------
col = [0.92 0.41 0.20; 0.38 0.31 0.84; 0.16 0.47 0.84];           % A orange, B violet, C blue
anchors = [255 255 255; 207 225 247; 127 176 234; 42 120 214; 18 63 122; 7 29 58]/255;
cmap = interp1(linspace(0, 1, 6), anchors, linspace(0, 1, 256));
lim = 11.5; sel = abs(x) <= lim;

figure(1); clf;                               % refractive-index maps
for q = 1:3
    subplot(1, 3, q); imagesc(x(sel), x(sel), sqrt(maps{q}(sel, sel))); axis image; set(gca, 'YDir', 'normal');
    caxis([1 1.45]); xlabel('x (\mum)'); ylabel('y (\mum)'); title(names{q});
end
colormap(flipud(gray(256)));
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 4]); print(gcf, 'bpm_fig1_index_maps.png', '-dpng', '-r200');

figure(2); clf;                               % intensity at three distances, in dB below the launched peak
for q = 1:3
    for j = 1:3
        subplot(3, 3, 3*(q - 1) + j);
        imagesc(x(sel), x(sel), 10*log10(res{q}.snap{j}(sel, sel) + 1e-12)); axis image; set(gca, 'YDir', 'normal');
        caxis([-30 0]); title(sprintf('%s,  z = %.0f \\mum', des(q), res{q}.snap_z(j)));
        if q == 3, xlabel('x (\mum)'); end
        if j == 1, ylabel('y (\mum)'); end
    end
end
colormap(cmap);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 9.5 9]); print(gcf, 'bpm_fig2_intensity.png', '-dpng', '-r200');

figure(3); clf;                               % power in the liquid core against distance
for q = 1:3
    semilogy(res{q}.z, res{q}.Pcore, '-', 'Color', col(q, :), 'LineWidth', 1.8); hold on;
end
semilogy(res{1}.z, exp(polyval(pA, res{1}.z)), 'k--', 'LineWidth', 1.0);
semilogy(res{2}.z, res{2}.Ptot, ':', 'Color', col(2, :), 'LineWidth', 1.4);
hold off; grid on; axis([0 Z 1e-3 2]);
xlabel('distance z (\mum)'); ylabel('power in the liquid core / launched power');
legend('A: liquid core in solid silica', 'B: conventional cladding', 'C: this work', sprintf('fit for A: %.0f dB/mm', aA*dB), ...
       'B: all power left in the window', 'Location', 'southwest');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 7 4.6]); print(gcf, 'bpm_fig3_core_power.png', '-dpng', '-r200');

figure(4); clf;                               % phase difference between the two glucose levels
plot(zg, dphi, '-', 'Color', col(3, :), 'LineWidth', 1.8); hold on;
plot(zg, polyval(pg, zg), 'k--', 'LineWidth', 1.0); hold off; grid on;
xlabel('distance z (\mum)'); ylabel('phase difference (rad)');
title(sprintf('design C, %d against %d mg/dL glucose', Cg(2), Cg(1)));
legend('beam propagation', sprintf('straight line, %.2e rad/\\mum', pg(1)), 'Location', 'northwest');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 4.2]); print(gcf, 'bpm_fig4_glucose_phase.png', '-dpng', '-r200');
