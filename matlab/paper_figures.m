% PAPER_FIGURES  Draw Figs. 5-12 of the paper in MATLAB from the computed data, as editable figures.
%
% Reads ../out/figure_data.mat (written by code/export_figure_data.py from the result files in out/) and draws
%   Fig. 5  stochastic verification          Fig. 9   beam propagation in three liquid-core fibers
%   Fig. 6  fundamental mode                  Fig. 10  side view of the beam propagation
%   Fig. 7  dispersion and effective area     Fig. 11  sensing
%   Fig. 8  propagation with loss             Fig. 12  glucose data set and accuracy test
% Each figure is saved in figures_paper/ as .fig (open it in MATLAB and edit colours, labels, limits, fonts ...)
% and as .png (300 dpi).  Figs. 2-4 are drawn by m1_phase_portraits.m and m2_exact_solutions.m.
% Runs in a few seconds; base MATLAB only (R2016b or later), also GNU Octave (there without the .fig files).
clear; close all;
F = load(fullfile('..', 'out', 'figure_data.mat'));
outdir = 'figures_paper'; if ~exist(outdir, 'dir'), mkdir(outdir); end
isOct = exist('OCTAVE_VERSION', 'builtin') > 0;

% colours of the paper
blue = [42 120 214]/255; orange = [235 104 52]/255; aqua = [27 175 122]/255; violet = [98 80 214]/255;
gray = [138 137 132]/255; ink = [11 11 11]/255; gridc = [217 216 212]/255; liq = [127 176 234]/255;
anchors = [255 255 255; 207 225 247; 127 176 234; 42 120 214; 18 63 122; 7 29 58]/255;
seq = interp1(linspace(0, 1, 6), anchors, linspace(0, 1, 256));         % white -> dark blue (intensity maps)
set(0, 'DefaultAxesFontSize', 10, 'DefaultLineLineWidth', 1.4, 'DefaultAxesBox', 'on', 'DefaultAxesTickDir', 'in');
cols3 = {blue, orange, aqua}; mk3 = {'o', 's', '^'};
% ---- how each figure is stored ----
figs = {}; names = {};

%% ---------------- Fig. 5: stochastic verification ----------------
d = F.f5; h = figure('Color', 'w', 'Position', [40 40 1000 680], 'Name', 'Fig. 5');
subplot(2, 2, 1);
plot(d.t, d.re_det, 'Color', gray, 'LineWidth', 1.0); hold on;
plot(d.t, d.re_num, 'Color', blue, 'LineWidth', 1.2);
plot(d.t(1:12:end), d.re_ex(1:12:end), 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ink);
plot(d.t, d.I_num, '--', 'Color', orange); hold off;
xlim([0 10]); ylim([-1.25 1.9]); xlabel('t'); ylabel('field at the soliton center');
legend({'Re L(0,t), \rho = 0', 'Re L(0,t), numerical, \rho = 0.3', 'Re L(0,t), exact solution', '|L(0,t)|^2, numerical'}, ...
       'Location', 'north', 'FontSize', 8, 'Box', 'off');
title('(a)  one noise realization');
subplot(2, 2, 2);
loglog(d.bright(:, 1), d.bright(:, 2), 'o-', 'Color', blue, 'MarkerFaceColor', 'w'); hold on;
loglog(d.dark(:, 1), d.dark(:, 2), 's-', 'Color', orange, 'MarkerFaceColor', 'w');
dd = [6e-4 1.2e-2]; loglog(dd, 0.9*dd, ':', 'Color', gray, 'LineWidth', 0.8); hold off; grid on;
xlabel('time step \Deltat'); ylabel('max_x |L_{num} - L_{exact}| at t = 10');
legend({sprintf('bright soliton (slope %.2f)', d.bright_order), sprintf('dark soliton (slope %.2f)', d.dark_order), 'slope 1'}, ...
       'Location', 'northwest', 'Box', 'off');
title('(b)  pathwise convergence');
tt = linspace(0, 10, 200); k = 5:10:numel(d.te);
subplot(2, 2, 3); hold on; hm = zeros(1, 3); lbl = arrayfun(@(v) sprintf('\\rho = %.1f', v), d.rho, 'UniformOutput', false);
for q = 1:3
    plot(tt, exp(-d.rho(q)^2*tt/2), '-', 'Color', cols3{q}, 'LineWidth', 1.2);
    hm(q) = plot(d.te(k), d.meanfield(q, k), mk3{q}, 'Color', cols3{q}, 'MarkerFaceColor', 'w');
end
hold off; box on; xlim([0 10]); ylim([0.55 1.03]); xlabel('t'); ylabel('|\langle L(0,t)\rangle| / M(0)');
legend(hm, lbl, 'Location', 'southwest', 'Box', 'off'); title('(c)  mean field: lines exp(-\rho^2 t/2)');
subplot(2, 2, 4); hold on; hm = zeros(1, 3);
for q = 1:3
    plot(tt, d.rho(q)^2*tt, '-', 'Color', cols3{q}, 'LineWidth', 1.2);
    hm(q) = plot(d.te(k), d.phvar(q, k), mk3{q}, 'Color', cols3{q}, 'MarkerFaceColor', 'w');
end
hold off; box on; xlim([0 10]); ylim([0 1]); xlabel('t'); ylabel('phase variance (rad^2)');
legend(hm, lbl, 'Location', 'northwest', 'Box', 'off'); title('(d)  phase diffusion: lines \rho^2 t');
figs{end+1} = h; names{end+1} = 'Fig05_stochastic';

%% ---------------- Fig. 6: fundamental mode ----------------
d = F.f6; h = figure('Color', 'w', 'Position', [60 60 1000 430], 'Name', 'Fig. 6');
lim = 5.2; s = d.xc < lim + 0.3; xq = d.xc(s); Q = d.Sz(s, s);
X = [-fliplr(xq), xq]; M = [rot90(Q, 2), flipud(Q); fliplr(Q), Q];        % mirror the computed quadrant
subplot(1, 2, 1); imagesc(X, X, M'); set(gca, 'YDir', 'normal'); axis image; xlim([-lim lim]); ylim([-lim lim]);
colormap(gca, seq); caxis([0 1]); colorbar; xlabel('x (\mum)'); ylabel('y (\mum)');
th = linspace(0, 2*pi, 300); hold on; plot(2.45*cos(th), 2.45*sin(th), 'Color', gray, 'LineWidth', 0.5);
plot(2.53*cos(th), 2.53*sin(th), 'Color', gray, 'LineWidth', 0.5); hold off;
title('(a)  S_z of the fundamental mode (normalized)');
subplot(1, 2, 2);
semilogy(d.xc, abs(d.Sz(:, 1)) + 1e-30, 'Color', blue); hold on;
semilogy(d.yc, abs(d.Sz(1, :)) + 1e-30, '--', 'Color', orange);
yl = [1e-22 5]; patch([0 2.45 2.45 0], [yl(1) yl(1) yl(2) yl(2)], liq, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot([2.53 2.53], yl, 'Color', gray, 'LineWidth', 0.5, 'HandleVisibility', 'off'); hold off;
xlim([0 11.9]); ylim(yl); xlabel('distance from the fiber axis (\mum)'); ylabel('|S_z| / S_{z,max}');
legend({'along x (\Gamma-K direction)', 'along y (\Gamma-M direction)'}, 'Location', 'northeast', 'Box', 'off');
title('(b)  radial decay');
figs{end+1} = h; names{end+1} = 'Fig06_mode';

%% ---------------- Fig. 7: dispersion and effective area ----------------
d = F.f7; h = figure('Color', 'w', 'Position', [80 80 1150 380], 'Name', 'Fig. 7');
subplot(1, 3, 1);
plot(d.lam, d.n_a, ':', 'Color', gray); hold on; plot(d.lam, d.neff, 'Color', blue); plot(d.lam, d.ng, '--', 'Color', orange); hold off;
xlim([0.86 1.34]); xlabel('wavelength (\mum)'); ylabel('index'); legend({'n_a (analyte)', 'n_{eff}', 'n_g'}, 'Location', 'east', 'Box', 'off'); title('(a)');
subplot(1, 3, 2);
plot([0.86 1.34], [0 0], 'Color', gridc, 'LineWidth', 0.6); hold on; plot(d.lam, d.beta2, 'Color', blue);
plot(d.ZDW, 0, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ink); plot(1.064, d.beta2_design, 's', 'Color', orange, 'MarkerFaceColor', orange); hold off;
text(d.ZDW + 0.03, 3, sprintf('ZDW %.3f \\mum', d.ZDW)); text(0.90, -58, 'design point'); text(0.90, -66, sprintf('%.1f ps^2/km', d.beta2_design));
xlim([0.86 1.34]); xlabel('wavelength (\mum)'); ylabel('\beta_2 (ps^2/km)'); title('(b)');
subplot(1, 3, 3);
plot(d.lam, d.Aeff, 'Color', blue); hold on; plot(1.064, d.Aeff_design, 's', 'Color', orange, 'MarkerFaceColor', orange); hold off;
xlim([0.86 1.34]); xlabel('wavelength (\mum)'); ylabel('A_{eff} (\mum^2)'); title('(c)');
figs{end+1} = h; names{end+1} = 'Fig07_dispersion';

%% ---------------- Fig. 8: propagation with loss and third-order dispersion ----------------
d = F.f8; h = figure('Color', 'w', 'Position', [100 100 1150 380], 'Name', 'Fig. 8');
T = {d.T30, d.T20}; lab = {'(a)  T_0 = 30 fs', '(b)  T_0 = 20 fs'};
for q = 1:2
    subplot(1, 3, q); r = T{q};
    plot([5 5], [0.9 4.2], 'Color', gridc, 'LineWidth', 0.6, 'HandleVisibility', 'off'); hold on;
    plot(r.z_cm, r.linear, '--', 'Color', orange); plot(r.z_cm, r.soliton, 'Color', blue); plot(r.z_cm, r.ideal, ':', 'Color', gray); hold off;
    xlim([0 8]); ylim([0.9 4.2]); xlabel('distance z (cm)'); ylabel('FWHM / input FWHM'); title(lab{q});
    if q == 1, legend({'linear (dispersion + loss)', 'soliton + loss + TOD', 'ideal soliton'}, 'Location', 'northwest', 'Box', 'off'); end
end
subplot(1, 3, 3); s = d.snap;
plot(s.T_fs, sech(s.T_fs/30).^2, ':', 'Color', gray); hold on; plot(s.T_fs, s.soliton, 'Color', blue); plot(s.T_fs, s.linear, '--', 'Color', orange); hold off;
xlim([-200 200]); ylim([0 1.05]); xlabel('retarded time T (fs)'); ylabel('power / input peak power');
legend({'input', 'soliton, z = 5 cm', 'linear, z = 5 cm'}, 'Location', 'northeast', 'Box', 'off'); title('(c)  T_0 = 30 fs');
figs{end+1} = h; names{end+1} = 'Fig08_propagation';

%% ---------------- Fig. 9: beam propagation in three liquid-core fibers ----------------
d = F.f9; h = figure('Color', 'w', 'Position', [120 40 1250 760], 'Name', 'Fig. 9');
idxmap = interp1([1 1 + 0.34 1.45], [1 1 1; liq; 0.79 0.78 0.76], linspace(1, 1.45, 256));   % air white, liquid blue, silica gray
dnames = {'A  liquid core in solid silica', 'B  conventional cladding, d/\Lambda = 0.5', 'C  this work, d/\Lambda = 0.952'};
w = 0.105; g = 0.008; x0 = 0.05;
for q = 1:3
    y0 = 0.68 - (q - 1)*0.30;
    for j = 1:4
        ax = axes('Position', [x0 + (j - 1)*(w + g), y0, w, w*1250/760]);
        if j == 1
            imagesc(d.x, d.x, squeeze(d.index(q, :, :))); colormap(ax, idxmap); caxis([1 1.45]);
            text(-10.5, -10, 'index', 'FontSize', 8);
        else
            imagesc(d.x, d.x, squeeze(d.intensity_dB(q, j - 1, :, :))); colormap(ax, seq); caxis([-30 0]);
            title(sprintf('z = %d \\mum', d.snap_z(j - 1)), 'FontSize', 9, 'FontWeight', 'normal');
        end
        set(ax, 'YDir', 'normal', 'XTick', [-8 0 8], 'YTick', [-8 0 8], 'FontSize', 8); axis(ax, 'image');
        if q < 3, set(ax, 'XTickLabel', []); else, xlabel('x (\mum)'); end
        if j > 1, set(ax, 'YTickLabel', []); else, ylabel('y (\mum)'); end
    end
    annotation('textbox', [x0, y0 + w*1250/760 + 0.012, 0.4, 0.03], 'String', dnames{q}, 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 9);
end
annotation('textbox', [0.01 0.95 0.05 0.04], 'String', '(a)', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 11);
cb = colorbar(ax, 'southoutside'); set(cb, 'Position', [x0 + w + g, 0.045, 3*w + 2*g, 0.015]); xlabel(cb, 'intensity (dB below the launched peak)');
ax = axes('Position', [0.60 0.56 0.37 0.37]); cc = {orange, violet, blue}; hm = zeros(1, 4);
for q = [3 2 1], hm(4 - q) = semilogy(d.z, d.Pcore(q, :), 'Color', cc{q}); hold on; end
hm(4) = semilogy(d.z, exp(-d.decay_A_exact*d.z), 'k--', 'LineWidth', 1.0); hold off;
xlim([0 200]); ylim([1e-3 2]); xlabel('distance z (\mum)'); ylabel('power in the liquid core / launched power');
legend(hm, {'C: this work', 'B: conventional cladding', 'A: liquid core in solid silica', 'exact leaky mode of A'}, 'Location', 'southwest', 'Box', 'off');
title('(b)');
ax = axes('Position', [0.60 0.08 0.37 0.36]);
plot(d.gz, d.gdphi, 'Color', blue); xlim([0 100]); xlabel('distance z (\mum)'); ylabel('phase difference (rad)');
text(5, 0.175, 'design C: 300 against 80 mg/dL glucose', 'FontSize', 9);
text(5, 0.160, sprintf('slope %.2f mrad/\\mum = %.1f rad over 5 cm', 1e3*d.gslope, 5e4*d.gslope), 'FontSize', 9);
title('(c)');
figs{end+1} = h; names{end+1} = 'Fig09_beam_propagation';

%% ---------------- Fig. 10: side view of the beam propagation ----------------
d = F.f10; h = figure('Color', 'w', 'Position', [140 140 1250 400], 'Name', 'Fig. 10');
ttl = {'(a)  A: liquid core in solid silica', '(b)  B: conventional cladding, d/\Lambda = 0.5', '(c)  C: this work, d/\Lambda = 0.952'};
for q = 1:3
    ax = subplot(1, 3, q);
    imagesc(d.zs, d.xs, d.side_dB(:, :, q)); set(ax, 'YDir', 'normal'); colormap(ax, seq); caxis([-30 0]);
    hold on; plot([0 200], [d.Rc d.Rc], 'k:', [0 200], -[d.Rc d.Rc], 'k:', 'LineWidth', 0.8); hold off;
    xlabel('distance z (\mum)'); if q == 1, ylabel('x (\mum), on the line y = 0'); end
    title(ttl{q});
end
cb = colorbar; ylabel(cb, 'intensity (dB below launched peak)');
figs{end+1} = h; names{end+1} = 'Fig10_side_view';

%% ---------------- Fig. 11: sensing ----------------
d = F.f11; h = figure('Color', 'w', 'Position', [160 40 1000 680], 'Name', 'Fig. 11');
subplot(2, 2, 1); p = polyfit(d.n_a, d.neff, 1); xx = linspace(min(d.n_a), max(d.n_a), 50);
plot(xx, polyval(p, xx), 'Color', blue, 'LineWidth', 1.2); hold on;
plot(d.n_a, d.neff, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ink); hold off;
xlabel('analyte index n_a'); ylabel('n_{eff}'); legend({sprintf('linear fit, slope %.4f', p(1)), 'mode solver'}, 'Location', 'northwest', 'Box', 'off');
title('(a)  index sensitivity');
Cg = linspace(0, 300, 301); Ls = [2 5 10]; ls = {'-.', '-', '--'}; lc = {aqua, blue, orange};
Llbl = arrayfun(@(v) sprintf('L_f = %g cm', v), Ls, 'UniformOutput', false);
subplot(2, 2, 2); hold on;
for q = 1:3, plot(Cg, d.dphi_per_mgdl*Cg*Ls(q)/5, ls{q}, 'Color', lc{q}); end
hold off; box on; xlim([0 300]); ylim([0 260]); xlabel('glucose concentration change (mg/dL)'); ylabel('phase shift \Delta\Phi (rad)');
legend(Llbl, 'Location', 'northwest', 'Box', 'off'); title('(b)  phase read-out');
subplot(2, 2, 3); hold on;
for q = 1:3, plot(Cg, d.dtau_fs_per_mgdl*Cg*Ls(q)/5, ls{q}, 'Color', lc{q}); end
hold off; box on; xlim([0 300]); ylim([0 145]); xlabel('glucose concentration change (mg/dL)'); ylabel('group-delay shift \Delta\tau (fs)');
legend(Llbl, 'Location', 'northwest', 'Box', 'off'); title('(c)  arrival-time read-out');
subplot(2, 2, 4); vr = logspace(-3, log10(0.5), 100); k0 = 2*pi/d.lam0_m;
for q = 1:3
    L = Ls(q)/100; lod = 3*vr*sqrt(L/d.LD_m)/(k0*L*d.Sn*d.dndC);
    hm(q) = loglog(vr, lod, ls{q}, 'Color', lc{q}); hold on;
end
loglog([1e-3 0.5], [15 15], ':', 'Color', gray, 'LineWidth', 0.8); hold off; grid on;
text(1.15e-3, 19, 'ISO 15197:2013 limit (15 mg/dL)', 'FontSize', 8, 'Color', gray);
xlim([1e-3 0.5]); ylim([5e-3 40]); xlabel('noise strength \rho'); ylabel('3\sigma detection limit (mg/dL)');
legend(hm, Llbl, 'Location', 'southeast', 'Box', 'off'); title('(d)  noise-limited detection limit');
figs{end+1} = h; names{end+1} = 'Fig11_sensing';

%% ---------------- Fig. 12: glucose data set and virtual accuracy test ----------------
d = F.f12; h = figure('Color', 'w', 'Position', [180 180 1200 420], 'Name', 'Fig. 12');
subplot(1, 3, 1); xx = 1:numel(d.levels);
b = bar(xx, [d.counts(:) - d.counts_diabetes(:), d.counts_diabetes(:)]/1000, 0.72, 'stacked', 'EdgeColor', 'w');
set(b(1), 'FaceColor', blue); set(b(2), 'FaceColor', orange);
set(gca, 'XTick', xx, 'XTickLabel', arrayfun(@(v) sprintf('%d', v), d.levels, 'UniformOutput', false), 'FontSize', 8);
if ~isOct, xtickangle(90); end
ylim([0 10.5]); xlabel('blood glucose level in the data set (mg/dL)'); ylabel('records (thousands)');
legend({'no diabetes label', 'diabetes label'}, 'Location', 'northeast', 'Box', 'off'); title('(a)  test concentrations');
subplot(1, 3, 2); i0 = find(d.levels == 100); dn = d.neff - d.neff(i0); p = polyfit(d.levels, dn, 1); cc = linspace(70, 310, 50);
plot(cc, 1e4*polyval(p, cc), 'Color', blue, 'LineWidth', 1.2); hold on; plot(d.levels, 1e4*dn, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ink); hold off;
xlim([60 320]); xlabel('glucose concentration (mg/dL)'); ylabel('n_{eff} - n_{eff}(100 mg/dL)  (10^{-4})');
legend({'linear fit', 'mode solver'}, 'Location', 'northwest', 'Box', 'off'); title('(b)  computed response');
subplot(1, 3, 3);
plot(d.sigT, 100*d.within_rho01, 'Color', blue); hold on; plot(d.sigT, 100*d.within_rho10, '--', 'Color', orange);
plot([0 0.4], [95 95], ':', 'Color', gray, 'LineWidth', 0.8); hold off;
text(0.39, 97, 'ISO 15197:2013 requirement (95 %)', 'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', gray);
xlim([0 0.4]); ylim([60 101.5]); xlabel('temperature uncertainty \sigma_\vartheta (K)'); ylabel('records within ISO limits (%)');
legend({'\rho = 0.1', '\rho = 1.0'}, 'Location', 'southwest', 'Box', 'off'); title('(c)  virtual accuracy test');
figs{end+1} = h; names{end+1} = 'Fig12_dataset';

%% ---------------- save ----------------
for k = 1:numel(figs)
    f = fullfile(outdir, names{k});
    if ~isOct, savefig(figs{k}, [f '.fig']); end
    set(figs{k}, 'PaperPositionMode', 'auto'); print(figs{k}, [f '.png'], '-dpng', '-r300');
end
fprintf('saved %d figures in %s (.fig for editing in MATLAB, .png at 300 dpi)\n', numel(figs), fullfile(pwd, outdir));
