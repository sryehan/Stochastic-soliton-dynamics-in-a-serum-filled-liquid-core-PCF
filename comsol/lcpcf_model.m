function model = lcpcf_model(varargin)
%LCPCF_MODEL  COMSOL model of the serum-filled liquid-core photonic crystal fiber of the paper
%
%   S. Ryehan et al., "Stochastic soliton dynamics in a serum-filled liquid-core photonic crystal fiber
%   for refractometric glucose sensing", Section 6.
%
%   Runs through LiveLink for MATLAB:
%       >> model = lcpcf_model;                     % fundamental mode at 1064 nm
%       >> model = lcpcf_model('n_a', 1.34392);     % a second index, for the sensitivity d n_eff / d n_a
%       >> model = lcpcf_model('lam0', 1.10);       % another wavelength, for the dispersion
%       >> model = lcpcf_model('nrings', 2);        % fewer rings, for the confinement loss
%   Start MATLAB from the COMSOL "with MATLAB" shortcut, or run  comsol mphserver matlab  first.
%
%   *** THIS SCRIPT HAS NOT BEEN EXECUTED BY ITS AUTHOR: no COMSOL licence was available. ***
%   The structure is that of a COMSOL Model M-file, but a few property names differ between COMSOL
%   versions.  If a line fails, build that one step in the COMSOL Desktop (COMSOL_guide.md lists every
%   step with its numbers) and use File > Save As > Model File for MATLAB to see the syntax of your
%   version.  The expected results are printed at the end and listed in COMSOL_guide.md.

import com.comsol.model.*
import com.comsol.model.util.*

% ------------------------------------------------------------------ parameters
p = struct('lam0', 1.064, ...     % wavelength (um)
           'Lambda', 1.66, ...    % pitch (um)
           'dh', 1.58, ...        % air-hole diameter (um)
           'Rc', 2.45, ...        % radius of the liquid core (um)
           'tw', 0.08, ...        % thickness of the silica core wall (um)
           'nrings', 5, ...       % number of rings of holes (rings 2 .. nrings+1 of the lattice)
           'n_a', [], ...         % analyte index; default = water(lam0) + 0.016
           'dpml', 2.5, ...       % PML thickness (um)
           'neigs', 10, ...       % number of modes to compute
           'hwall', 0.03, ...     % mesh size in the core wall (um)
           'hcore', 0.15, ...     % mesh size in the core and the inner cladding (um)
           'hout', 0.50);         % mesh size in the outer cladding (um)
for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end
n_sio2 = sellmeier_silica(p.lam0);
n_w    = kedenburg_water(p.lam0);
if isempty(p.n_a), p.n_a = n_w + 0.016; end            % serum = water + 0.016 (Section 6.1)
Rjacket = p.nrings*p.Lambda + p.dh/2 + 1.0;            % outer radius of the silica region
fprintf('lam0 = %.4f um: n_silica = %.5f, n_water = %.5f, n_a = %.5f\n', p.lam0, n_sio2, n_w, p.n_a);

model = ModelUtil.create('LCPCF');
model.label('liquid-core PCF, serum-filled');
model.param.set('lam0', sprintf('%.6f[um]', p.lam0));
model.param.set('f0', 'c_const/lam0');
model.param.set('n_sio2', sprintf('%.6f', n_sio2));
model.param.set('n_a', sprintf('%.6f', p.n_a));
comp = model.component.create('comp1', true);
geom = comp.geom.create('geom1', 2);
geom.lengthUnit('um');

% ------------------------------------------------------------------ geometry
% jacket, with the PML as an outer layer of the same circle
cj = geom.create('cjack', 'Circle');
cj.set('r', Rjacket + p.dpml);
cj.set('layername', {'pml'});
cj.setIndex('layer', p.dpml, 0);
% air holes: the triangular-lattice sites of rings 2 .. nrings+1 (ring 0 and ring 1 make room for the core)
xy = lattice_sites(p.Lambda, p.nrings);
hole_tags = cell(1, size(xy, 1));
for q = 1:size(xy, 1)
    tag = sprintf('h%03d', q); hole_tags{q} = tag;
    h = geom.create(tag, 'Circle');
    h.set('r', p.dh/2);
    h.set('pos', [xy(q,1) xy(q,2)]);
end
fprintf('%d air holes created (rings 2 to %d)\n', size(xy, 1), p.nrings + 1);
% the core wall truncates the six holes at sqrt(3)*Lambda; the six at 2*Lambda only touch it
cwall = geom.create('cwall', 'Circle'); cwall.set('r', p.Rc + p.tw);
dif = geom.create('trim', 'Difference');
dif.selection('input').set(hole_tags);
dif.selection('input2').set({'cwall'});
ccore = geom.create('ccore', 'Circle'); ccore.set('r', p.Rc);
geom.run;

% ------------------------------------------------------------------ domain selections (by coordinates)
d_core = mphselectcoords(model, 'geom1', [0; 0], 'domain', 'radius', 1e-4);
d_wall = mphselectcoords(model, 'geom1', [p.Rc + p.tw/2; 0], 'domain', 'radius', 1e-4);
d_pml  = mphselectcoords(model, 'geom1', [Rjacket + p.dpml/2; 0], 'domain', 'radius', 1e-4);
d_air  = [];
for q = 1:size(xy, 1)                                   % every hole centre is still inside its own domain
    d_air = [d_air, mphselectcoords(model, 'geom1', [xy(q,1); xy(q,2)], 'domain', 'radius', 1e-4)]; %#ok<AGROW>
end
d_air = unique(d_air);
fprintf('domains: core %d, wall %d, air %d, PML %d (of %d in all)\n', numel(d_core), numel(d_wall), numel(d_air), numel(d_pml), geom.getNDomains);
assert(numel(d_core) == 1 && numel(d_wall) == 1 && numel(d_air) == size(xy, 1), 'a selection did not find the expected domains - check the geometry');
sel = @(tag, dom) make_sel(comp, tag, dom);
sel('sel_core', d_core); sel('sel_air', d_air); sel('sel_pml', d_pml);

% ------------------------------------------------------------------ materials (assigned in this order)
add_material(comp, 'silica', 'n_sio2', []);             % everything, then overridden below
add_material(comp, 'air',    '1',      'sel_air');
add_material(comp, 'serum',  'n_a',    'sel_core');

% ------------------------------------------------------------------ PML, physics
pml = comp.coordSystem.create('pml1', 'PML');
pml.set('ScalingType', 'Cylindrical');
pml.selection.named('sel_pml');
ewfd = comp.physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
ewfd.prop('components').set('components', 'threeComponent');
ewfd.feature('wee1').set('DisplacementFieldModel', 'RefractiveIndex');

% ------------------------------------------------------------------ mesh
mesh = comp.mesh.create('mesh1');
ft = mesh.create('ftri1', 'FreeTri');
s0 = mesh.feature('size'); s0.set('hauto', 3); s0.set('custom', 'on'); s0.set('hmax', p.hout); s0.set('hmaxactive', 'on');
s1 = ft.create('sz_core', 'Size'); s1.selection.geom('geom1', 2); s1.selection.named('sel_core');
s1.set('custom', 'on'); s1.set('hmax', p.hcore); s1.set('hmaxactive', 'on');
s2 = ft.create('sz_wall', 'Size'); s2.selection.geom('geom1', 2); s2.selection.set(d_wall);
s2.set('custom', 'on'); s2.set('hmax', p.hwall); s2.set('hmaxactive', 'on');
mesh.run;
fprintf('mesh: %d elements\n', mphmesh_stats(model));

% ------------------------------------------------------------------ study: mode analysis
std = model.study.create('std1');
ma = std.create('mode', 'ModeAnalysis');
ma.set('modeFreq', 'f0');
ma.set('neigsactive', true); ma.set('neigs', p.neigs);
ma.set('shiftactive', true); ma.set('shift', sprintf('%.5f', p.n_a - 0.008));   % search for modes around this index
std.run;

% ------------------------------------------------------------------ results
neff = mphglobal(model, 'ewfd.neff', 'dataset', 'dset1', 'solnum', 'all');
i_all = model.result.numerical.create('i_all', 'IntSurface');
i_all.selection.all; i_all.set('expr', {'ewfd.Poavz', 'ewfd.Poavz^2'});
i_cor = model.result.numerical.create('i_cor', 'IntSurface');
i_cor.selection.named('sel_core'); i_cor.set('expr', {'ewfd.Poavz'});
fprintf('\n  mode   n_eff (real)     Im(n_eff)    core power     A_eff (um^2)\n');
best = struct('frac', -1);
for q = 1:numel(neff)
    i_all.set('solnum', q); i_cor.set('solnum', q);
    a = i_all.getReal(); c = i_cor.getReal();
    frac = c(1)/a(1); Aeff = a(1)^2/a(2)*1e12;          % m^2 -> um^2
    fprintf('   %2d   %12.7f  %11.3e     %7.3f %%     %8.3f\n', q, real(neff(q)), imag(neff(q)), 100*frac, Aeff);
    if frac > best.frac && real(neff(q)) < p.n_a
        best = struct('frac', frac, 'neff', neff(q), 'Aeff', Aeff, 'mode', q);
    end
end
loss_dB_m = 8.686*2*pi/(p.lam0*1e-6)*imag(best.neff);
fprintf('\nfundamental mode (largest fraction of the power in the core): mode %d\n', best.mode);
fprintf('  n_eff          = %.7f\n  core power     = %.3f %%\n  A_eff          = %.3f um^2\n  confinement loss = %.3e dB/m\n', ...
        real(best.neff), 100*best.frac, best.Aeff, loss_dB_m);
if abs(p.lam0 - 1.064) < 1e-9 && abs(p.n_a - 1.33992) < 1e-4 && p.nrings == 5
    fprintf(['\ncompare with the finite-difference solver of the paper (Table 6):\n' ...
             '  n_eff      1.332018   (difference %.1e)\n' ...
             '  core power 99.19 %%     (difference %.2f points)\n' ...
             '  A_eff      11.78 um^2  (difference %.1f %%)\n' ...
             'Agreement to about 1e-4 in n_eff, 0.2 points in the power fraction and 2 %% in A_eff\n' ...
             'counts as a confirmation; the two methods use different meshes and different formulations.\n'], ...
            abs(real(best.neff) - 1.332018), abs(100*best.frac - 99.19), 100*abs(best.Aeff/11.78 - 1));
end
fprintf(['\nnext steps (see COMSOL_guide.md):\n' ...
         '  sensitivity  S_n  : run again with ''n_a'', %.5f and divide the change of n_eff by 0.004 (paper: 0.996)\n' ...
         '  dispersion   beta2: run for lam0 = 0.90 : 0.02 : 1.30 and use beta2 = lam^3/(2*pi*c^2) * d2n_eff/dlam^2\n' ...
         '                      (paper: -26.6 ps^2/km at 1.064 um, zero-dispersion wavelength 0.927 um)\n' ...
         '  confinement  loss : run with ''nrings'', 1 and 2 (paper: 10.5 dB/m and 2.6e-5 dB/m)\n'], p.n_a + 0.004);
end

% =========================================================================================== helpers
function xy = lattice_sites(Lambda, nrings)
%LATTICE_SITES  Centres of the air holes: triangular-lattice sites of rings 2 .. nrings+1.
xy = [];
for i = -(nrings+1):(nrings+1)
    for j = -(nrings+1):(nrings+1)
        r = max([abs(i), abs(j), abs(i + j)]);
        if r >= 2 && r <= nrings + 1
            xy(end+1, :) = [(i + j/2)*Lambda, j*sqrt(3)/2*Lambda]; %#ok<AGROW>
        end
    end
end
end

function make_sel(comp, tag, dom)
s = comp.selection.create(tag, 'Explicit');
s.geom(2); s.set(dom); s.label(tag);
end

function add_material(comp, tag, n_expr, selname)
m = comp.material.create(tag, 'Common');
pg = m.propertyGroup.create('RefractiveIndex', 'Refractive index');
pg.set('n', {n_expr});                 % if your version wants the full tensor, use the 9-component form:
pg.set('ki', {'0'});                   % {n n n} on the diagonal and '0' elsewhere
if isempty(selname), m.selection.all; else, m.selection.named(selname); end
end

function n = sellmeier_silica(l)
%SELLMEIER_SILICA  Malitson, J. Opt. Soc. Am. 55 (1965) 1205; l in um.
n = sqrt(1 + 0.6961663*l^2/(l^2 - 0.0684043^2) + 0.4079426*l^2/(l^2 - 0.1162414^2) + 0.8974794*l^2/(l^2 - 9.896161^2));
end

function n = kedenburg_water(l)
%KEDENBURG_WATER  Kedenburg et al., Opt. Mater. Express 2 (2012) 1588; l in um.
n = sqrt(1 + 0.75831*l^2/(l^2 - 0.01007) + 0.08495*l^2/(l^2 - 8.91377));
end

function n = mphmesh_stats(model)
try
    st = mphmeshstats(model, 'mesh1'); n = sum(st.numelem);
catch
    n = 0;
end
end
