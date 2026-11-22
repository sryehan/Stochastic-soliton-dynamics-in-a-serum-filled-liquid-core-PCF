function eps = pcf_index_map(design, x, n_si, n_a, s)
%PCF_INDEX_MAP  Relative permittivity n^2(x,y) of a liquid-core fiber cross-section on the grid x (um).
%   design 'A' : liquid core (radius 2.45 um) in solid silica, no holes (a liquid-filled capillary)
%   design 'B' : liquid core (radius 2.45 um) in a conventional cladding: pitch 4.46 um, hole diameter 2.23 um
%                (d/Lambda = 0.5), two rings of air holes
%   design 'C' : this work: pitch 1.66 um, hole diameter 1.58 um (d/Lambda = 0.952), core radius 2.45 um with an
%                80 nm silica wall, five rings of air holes (lattice rings 2-6, 120 holes)
%   n_si, n_a : refractive indices of silica and of the liquid in the core; air holes have n = 1.
%   Each pixel gets the average of n^2 over s x s sub-pixels (default s = 4), so that curved boundaries and the
%   thin core wall are represented on a grid that is coarser than the wall.
if nargin < 5, s = 4; end
N = numel(x); dx = x(2) - x(1);
xf = x(1) - dx/2 + ((1:N*s) - 0.5)*dx/s;                     % sub-pixel coordinates
ef = n_si^2*ones(N*s, N*s);
Rc = 2.45;
switch upper(design)
    case 'A'
        holes = zeros(0, 2); rh = 0; rwall = 0;
    case 'B'
        P = 4.46; rh = 2.23/2; rwall = 0; holes = lattice(P, 1, 2);
    case 'C'
        P = 1.66; rh = 1.58/2; rwall = Rc + 0.08; holes = lattice(P, 2, 6);
    otherwise
        error('unknown design');
end
for q = 1:size(holes, 1), ef = paint(ef, xf, holes(q,1), holes(q,2), rh, 1.0); end      % air holes
if rwall > 0, ef = paint(ef, xf, 0, 0, rwall, n_si^2); end                               % silica core wall
ef = paint(ef, xf, 0, 0, Rc, n_a^2);                                                     % liquid core
eps = reshape(mean(mean(reshape(ef, s, N, s, N), 1), 3), N, N);
end

function xy = lattice(P, r1, r2)
% triangular-lattice sites of the hexagonal rings r1 .. r2
xy = zeros(0, 2);
for i = -r2:r2
    for j = -r2:r2
        r = max([abs(i), abs(j), abs(i + j)]);
        if r >= r1 && r <= r2, xy(end+1, :) = [(i + j/2)*P, j*sqrt(3)/2*P]; end %#ok<AGROW>
    end
end
end

function ef = paint(ef, xf, x0, y0, r, val)
% set ef = val inside the circle of radius r centred at (x0, y0); rows are y, columns are x
ix = find(abs(xf - x0) <= r); iy = find(abs(xf - y0) <= r);
if isempty(ix) || isempty(iy), return; end
[X, Y] = meshgrid(xf(ix) - x0, xf(iy) - y0);
blk = ef(iy, ix); blk(X.^2 + Y.^2 <= r^2) = val; ef(iy, ix) = blk;
end
