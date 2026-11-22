% M7_PROPAGATION_MOVIE  Videos of the beam travelling along the three liquid-core fibers of m6_beam_propagation.m.
%
% The same Gaussian beam (1064 nm) is launched into
%   A  a liquid core in solid silica (capillary),
%   B  a liquid core in a conventional photonic-crystal cladding (pitch 4.46 um, d/Lambda = 0.5),
%   C  the fiber of this work (pitch 1.66 um, d/Lambda = 0.952),
% and propagated over 200 um with the same scalar beam propagation method as m6 (same grid, step, absorber).
% The three fibers are advanced TOGETHER, so the live window shows A, B and C at the same distance.
%
% Output, all in the folder  videos_m7/ :
%   A_solid_silica.mp4, B_conventional_cladding.mp4, C_this_work.mp4   one fiber per video (square frame)
%   ABC_comparison.mp4                                                 the three fibers side by side (wide frame)
%   side_view.png                                                      intensity on the line y = 0 against z
% Frames are at most 1920 x 1080 pixels (the limit of the MPEG-4 writer on Windows), 30 frames per second.
% Where MPEG-4 is not available (MATLAB on Linux), Motion-JPEG .avi files are written instead.
%
% Run time: about 6-10 min in MATLAB on a laptop (fast = true: about 60 %).  Needs pcf_index_map.m.
% Do not cover or minimise the figure windows while the videos are being written (getframe copies the screen).
clear; close all;
fast = false;                                 % true: dz = 0.05 um (about 1 % numerical loss per 100 um)
dzf  = 1;                                     % distance between video frames (um): 201 frames over 200 um
fps  = 30;                                    % frames per second
hold_end = 2;                                 % the last frame is shown for this many seconds
outdir = 'videos_m7';

% ---------------- parameters (as in m6_beam_propagation.m) ----------------
lam0 = 1.064; k0 = 2*pi/lam0;
n_si = 1.449631; n_w = 1.323917; n_a = n_w + 0.016; n0 = 1.332;
N = 640; dx = 0.05; x = (-N/2:N/2-1)*dx;
dz = 0.025; if fast, dz = 0.05; end
Z = 200; w0 = 1.9; Rc = 2.45; rabs = 11.0;
c0 = 299.792458; n_g = 1.3633;               % um/ps; group index of the mode of design C
des = 'ABC';
titles = {'A: liquid core in solid silica', 'B: conventional cladding, d/\Lambda = 0.5', 'C: this work, d/\Lambda = 0.952'};
vnames = {'A_solid_silica', 'B_conventional_cladding', 'C_this_work', 'ABC_comparison'};
footer = 'scalar beam propagation, \lambda = 1064 nm, serum-filled core; colour: intensity in dB below the launched peak';

anchors = [255 255 255; 207 225 247; 127 176 234; 42 120 214; 18 63 122; 7 29 58]/255;   % colours of Fig. 9
cm = interp1(linspace(0, 1, 6), anchors, linspace(0, 1, 256));
lim = 11.5; sel = find(abs(x) <= lim); xs = x(sel); ns = numel(sel);
j0 = find(abs(x) == min(abs(x)), 1);        % row of y = 0
every = round(dzf/dz); nst = round(Z/dz);
sideStep = 0.5; everyS = round(sideStep/dz); zs = (0:round(Z/sideStep))*sideStep;

% ---------------- operators ----------------
[X, Y] = meshgrid(x, x); R = hypot(X, Y);
kx = 2*pi/(N*dx)*[0:N/2-1, -N/2:-1]; [KX, KY] = meshgrid(kx, kx);
D  = exp(-1i*(KX.^2 + KY.^2)/(2*k0*n0)*dz);
Dh = exp(-1i*(KX.^2 + KY.^2)/(2*k0*n0)*dz/2);
sig = 1.0*max(0, (R - rabs)/(max(x) - rabs)).^2;
psi0 = exp(-R.^2/w0^2); I0 = max(abs(psi0(:)).^2);
maps = cell(1, 3); V = cell(1, 3); psi = cell(1, 3);
for q = 1:3
    maps{q} = pcf_index_map(des(q), x, n_si, n_a);
    V{q} = exp(1i*k0*(maps{q} - n0^2)/(2*n0)*dz - sig*dz);
    psi{q} = ifft2(Dh.*fft2(psi0));          % first half step of diffraction
end
side = zeros(ns, numel(zs), 3);

% ---------------- figures: wide (A | B | C) and square (one fiber) ----------------
dB0 = 10*log10(abs(psi0(sel, sel)).^2/I0 + 1e-12);
scr = get(0, 'ScreenSize'); fW = min(1600, 0.95*scr(3))/1600; fS = min(860, 0.85*scr(4))/860;   % fit the screen
figW = figure('Color', 'w', 'Units', 'pixels', 'Position', [20 40 round(1600*fW) round(640*fW)], 'Name', 'A | B | C', 'MenuBar', 'none');
for q = 1:3
    axW(q) = axes('Parent', figW, 'Position', [0.04 + (q - 1)*0.295, 0.15, 0.26, 0.68]); %#ok<SAGROW>
    hW(q) = imagesc(xs, xs, dB0, 'Parent', axW(q)); hold(axW(q), 'on'); %#ok<SAGROW>
    contour(axW(q), xs, xs, sqrt(maps{q}(sel, sel)), [1.2 1.39], 'LineColor', [0.6 0.6 0.6], 'LineWidth', 0.5);
    hold(axW(q), 'off'); axis(axW(q), 'image'); set(axW(q), 'YDir', 'normal', 'FontSize', 11); caxis(axW(q), [-30 0]);
    title(axW(q), titles{q}, 'FontSize', 13); xlabel(axW(q), 'x (\mum)'); if q == 1, ylabel(axW(q), 'y (\mum)'); end
end
for q = 1:3, colormap(axW(q), cm); end; cb = colorbar(axW(3)); set(cb, 'Position', [0.925 0.15 0.012 0.68]); ylabel(cb, 'dB');
tW = annotation(figW, 'textbox', [0 0.89 1 0.1], 'String', '', 'HorizontalAlignment', 'center', ...
                'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold');
annotation(figW, 'textbox', [0 0 1 0.05], 'String', footer, 'HorizontalAlignment', 'center', ...
           'EdgeColor', 'none', 'FontSize', 10, 'Color', [0.35 0.35 0.35]);

figS = figure('Color', 'w', 'Units', 'pixels', 'Position', [60 30 round(860*fS) round(860*fS)], 'Name', 'one fiber', 'MenuBar', 'none');
axS = axes('Parent', figS, 'Position', [0.1 0.12 0.72 0.72]);
hS = imagesc(xs, xs, dB0, 'Parent', axS); hold(axS, 'on');
[~, hc] = contour(axS, xs, xs, sqrt(maps{1}(sel, sel)), [1.2 1.39], 'LineColor', [0.6 0.6 0.6], 'LineWidth', 0.5);
hold(axS, 'off'); axis(axS, 'image'); set(axS, 'YDir', 'normal', 'FontSize', 12); caxis(axS, [-30 0]);
xlabel(axS, 'x (\mum)'); ylabel(axS, 'y (\mum)'); hT = title(axS, titles{1}, 'FontSize', 15);
colormap(axS, cm); cbS = colorbar(axS); set(cbS, 'Position', [0.85 0.12 0.025 0.72]); ylabel(cbS, 'intensity (dB)');
tS = annotation(figS, 'textbox', [0 0.9 1 0.08], 'String', '', 'HorizontalAlignment', 'center', ...
                'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold');
annotation(figS, 'textbox', [0.02 0.0 0.96 0.06], 'String', footer, 'HorizontalAlignment', 'center', ...
           'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.35 0.35 0.35]);
drawnow;

% frame sizes: shrink the windows until a captured frame fits into 1920 x 1080, then fix it (even numbers)
figs = [figS figW]; Hf = zeros(1, 2); Wf = zeros(1, 2);
for k = 1:2
    for attempt = 1:4
        im = getframe(figs(k)); im = im.cdata; s = min([1, 1920/size(im, 2), 1080/size(im, 1)]);
        if s >= 1, break; end
        p = get(figs(k), 'Position'); set(figs(k), 'Position', [p(1:2), floor(p(3:4)*s*0.98)]); drawnow;
    end
    Hf(k) = 2*floor(size(im, 1)/2); Wf(k) = 2*floor(size(im, 2)/2);
end
fprintf('video frames: %d x %d (single fiber), %d x %d (A | B | C)\n', Wf(1), Hf(1), Wf(2), Hf(2));

% ---------------- video writers ----------------
if ~exist(outdir, 'dir'), mkdir(outdir); end
prof = 'MPEG-4';
try
    pl = VideoWriter.getProfiles(); if ~any(strcmp({pl.Name}, prof)), prof = 'Motion JPEG AVI'; end
catch
    prof = 'Motion JPEG AVI';
end
vw = cell(1, 4);
for v = 1:4
    vw{v} = VideoWriter(fullfile(outdir, vnames{v}), prof);
    vw{v}.FrameRate = fps; vw{v}.Quality = 100;
    open(vw{v});
end
fprintf('writing %s videos to %s\n', prof, fullfile(pwd, outdir));

% ---------------- propagate A, B, C together and write the frames ----------------
Icur = cell(1, 3); for q = 1:3, Icur{q} = abs(psi0(sel, sel)).^2/I0; side(:, 1, q) = abs(psi0(j0, sel)').^2/I0; end
ks = 1; tic;
for it = 0:nst
    if it > 0
        rec = mod(it, every) == 0; recS = mod(it, everyS) == 0;
        for q = 1:3
            psi{q} = psi{q}.*V{q};
            if rec || recS
                f = ifft2(Dh.*fft2(psi{q})); I = abs(f).^2/I0;    % field at z = it*dz
                if recS, side(:, ks + 1, q) = I(j0, sel)'; end
                if rec, Icur{q} = I(sel, sel); end
            end
            psi{q} = ifft2(D.*fft2(psi{q}));
        end
        if recS, ks = ks + 1; end
        if ~rec, continue; end
    end
    z = it*dz; lbl = sprintf('z = %5.1f \\mum      travel time  t = %.3f ps', z, z*n_g/c0);
    nrep = 1; if it == nst, nrep = round(hold_end*fps); end
    % wide video: A | B | C
    for q = 1:3, set(hW(q), 'CData', 10*log10(Icur{q} + 1e-12)); end
    set(tW, 'String', lbl); drawnow;
    fr = getframe(figW); im = 255*ones(Hf(2), Wf(2), 3, 'uint8');
    h = min(Hf(2), size(fr.cdata, 1)); w = min(Wf(2), size(fr.cdata, 2)); im(1:h, 1:w, :) = fr.cdata(1:h, 1:w, :);
    for r = 1:nrep, writeVideo(vw{4}, im); end
    % single-fiber videos
    for q = 1:3
        set(hS, 'CData', 10*log10(Icur{q} + 1e-12));
        set(hc, 'ZData', sqrt(maps{q}(sel, sel))); set(hT, 'String', titles{q}); set(tS, 'String', lbl); drawnow;
        fr = getframe(figS); im = 255*ones(Hf(1), Wf(1), 3, 'uint8');
        h = min(Hf(1), size(fr.cdata, 1)); w = min(Wf(1), size(fr.cdata, 2)); im(1:h, 1:w, :) = fr.cdata(1:h, 1:w, :);
        for r = 1:nrep, writeVideo(vw{q}, im); end
    end
    if mod(round(z), 20) == 0, fprintf('  z = %3.0f um   (%.0f s)\n', z, toc); end
end
for v = 1:4, close(vw{v}); end
fprintf('done: 4 videos in %s (%.0f s)\n', fullfile(pwd, outdir), toc);

% ---------------- side view: intensity on y = 0 against distance ----------------
fig3 = figure('Color', 'w', 'Position', [100 100 1500 480]);
for q = 1:3
    subplot(1, 3, q);
    imagesc(zs, xs, 10*log10(side(:, :, q) + 1e-12)); set(gca, 'YDir', 'normal'); caxis([-30 0]);
    hold on; plot([0 Z], [Rc Rc], 'k:', [0 Z], -[Rc Rc], 'k:'); hold off;   % core boundary
    xlabel('distance z (\mum)'); if q == 1, ylabel('x (\mum), on the line y = 0'); end
    title(titles{q});
end
colormap(cm); colorbar('eastoutside');
set(fig3, 'PaperUnits', 'inches', 'PaperPosition', [0 0 15 4.8]);
print(fig3, fullfile(outdir, 'side_view.png'), '-dpng', '-r200');
fprintf('wrote %s\n', fullfile(outdir, 'side_view.png'));
