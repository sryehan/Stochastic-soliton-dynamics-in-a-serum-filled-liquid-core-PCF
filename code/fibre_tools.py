"""
Fibre tools: material data, exact step-index HE11 solver, wrapper around the full-vector
finite-difference mode solver (Fallahkhair, Li & Murphy, J. Lightwave Technol. 26, 1423 (2008))
as implemented in the open-source EMpy package, with sub-pixel permittivity averaging and an
optional complex-coordinate (PML) stretch for leaky modes.
"""
import numpy as np
from scipy.special import jv, kv, jvp, kvp
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')
from EMpy.modesolvers.FD import VFDModeSolver

c0 = 299792458.0

# ---------------- materials (lam in micrometres) ----------------
def n_silica(lam):
    """Malitson, J. Opt. Soc. Am. 55, 1205 (1965), fused silica, 20 C"""
    l2 = lam**2
    return np.sqrt(1 + 0.6961663*l2/(l2 - 0.0684043**2) + 0.4079426*l2/(l2 - 0.1162414**2)
                   + 0.8974794*l2/(l2 - 9.896161**2))

def n_water(lam):
    """Kedenburg et al., Opt. Mater. Express 2, 1588 (2012), distilled water, 20 C, 0.5-1.6 um"""
    l2 = lam**2
    return np.sqrt(1 + 0.75831*l2/(l2 - 0.01007) + 0.08495*l2/(l2 - 8.91377))

# ---------------- exact HE11 of a step-index fibre ----------------
def he11_exact(lam, a, n1, n2):
    k0 = 2*np.pi/lam
    V = k0*a*np.sqrt(n1**2 - n2**2)
    def F(ne):
        be = ne*k0
        u = a*np.sqrt(k0**2*n1**2 - be**2); w = a*np.sqrt(be**2 - k0**2*n2**2)
        J = jvp(1, u)/(u*jv(1, u)); K = kvp(1, w)/(w*kv(1, w))
        return (J + K)*(n1**2*J + n2**2*K) - ne**2*(1/u**2 + 1/w**2)**2
    # HE11 is the root with the largest n_eff: scan downward from n1
    ne = np.linspace(n1 - 1e-9, n2 + 1e-9, 20001)
    f = np.array([F(v) for v in ne])
    for i in range(len(ne) - 1):
        if np.isfinite(f[i]) and np.isfinite(f[i+1]) and f[i]*f[i+1] < 0 and abs(f[i]) < 1e6 and abs(f[i+1]) < 1e6:
            return brentq(F, ne[i+1], ne[i])
    raise RuntimeError('no root')

# ---------------- geometry helpers ----------------
def hex_centres(pitch, nrings):
    """centres of a triangular lattice, hexagonal outline, rings 1..nrings (centre excluded)"""
    pts = []
    for i in range(-nrings, nrings + 1):
        for j in range(-nrings, nrings + 1):
            if i == 0 and j == 0:
                continue
            if abs(i + j) > nrings:
                continue
            pts.append((pitch*(i + 0.5*j), pitch*np.sqrt(3)/2*j))
    return np.array(pts)

def averaged_eps(xc, yc, shapes, eps_bg, ss=8):
    """Arithmetic sub-pixel average of the permittivity on the cell-centre grid (xc, yc).
    shapes: list of (x0, y0, radius, eps).  Later shapes override earlier ones."""
    dxc = np.gradient(xc.real); dyc = np.gradient(yc.real)
    off = (np.arange(ss) + 0.5)/ss - 0.5
    X = (xc.real[:, None] + dxc[:, None]*off[None, :]).ravel()
    Y = (yc.real[:, None] + dyc[:, None]*off[None, :]).ravel()
    XX, YY = np.meshgrid(X, Y, indexing='ij')
    eps = np.full(XX.shape, eps_bg, dtype=complex)
    for (x0, y0, rad, e) in shapes:
        if (x0 + rad < X.min()) or (x0 - rad > X.max()) or (y0 + rad < Y.min()) or (y0 - rad > Y.max()):
            continue
        ix = np.where(np.abs(X - x0) <= rad)[0]; iy = np.where(np.abs(Y - y0) <= rad)[0]
        if len(ix) == 0 or len(iy) == 0:
            continue
        sub = (XX[np.ix_(ix, iy)] - x0)**2 + (YY[np.ix_(ix, iy)] - y0)**2 <= rad**2
        blk = eps[np.ix_(ix, iy)]; blk[sub] = e; eps[np.ix_(ix, iy)] = blk
    eps = eps.reshape(len(xc), ss, len(yc), ss).mean(axis=(1, 3))
    return eps if np.iscomplexobj(eps) and np.abs(eps.imag).max() > 0 else eps.real

def stretched(xr, x_pml, d_pml, smax):
    """complex coordinate x~ = x + i*int sigma, sigma = smax*((x-x_pml)/d_pml)^2 beyond x_pml"""
    xs = xr.astype(complex)
    m = xr > x_pml
    xs[m] = xr[m] + 1j*smax*d_pml*((xr[m] - x_pml)/d_pml)**3/3.0
    return xs

def solve_modes(lam, x, y, shapes, eps_bg, boundary='0A0A', neigs=6, n_target=None, ss=8, eps=None):
    """Shift-invert eigen-solve around n_target; returns list of dicts sorted by Re(n_eff) descending.
    Fields are returned on the node grid (x, y) for H and on the cell-centre grid for E."""
    from scipy.sparse import linalg
    from EMpy.modesolvers.FD import FDMode
    xc = (x[:-1] + x[1:])/2; yc = (y[:-1] + y[1:])/2
    if eps is None:
        eps = averaged_eps(xc, yc, shapes, eps_bg, ss)
    s = VFDModeSolver(lam, x, y, lambda a, b: eps, boundary)
    A = s.build_matrix().tocsc()
    k = 2*np.pi/lam
    vals, vecs = linalg.eigs(A, k=neigs, sigma=(n_target*k)**2, which='LM', tol=1e-10)
    neffs = np.sqrt(vals)/k
    idx = np.argsort(-neffs.real)
    neffs = neffs[idx]; vecs = vecs[:, idx]
    nx, ny = s.nx, s.ny
    Hxs = [vecs[:nx*ny, i].reshape(nx, ny) for i in range(neigs)]
    Hys = [vecs[nx*ny:, i].reshape(nx, ny) for i in range(neigs)]
    Hzs, Exs, Eys, Ezs = s.compute_other_fields(neffs, Hxs, Hys)
    modes = []
    for i in range(neigs):
        modes.append(dict(neff=neffs[i], Hx=Hxs[i], Hy=Hys[i], Hz=Hzs[i], Ex=Exs[i], Ey=Eys[i], Ez=Ezs[i]))
    return modes, eps, (xc, yc)

def cell_weights(x, y):
    """area elements on the cell-centre grid (real part of coordinates)"""
    return np.outer(np.diff(x.real), np.diff(y.real))

def centre(F):
    """average a node-grid field onto cell centres"""
    return 0.25*(F[:-1, :-1] + F[1:, :-1] + F[:-1, 1:] + F[1:, 1:])

def mode_integrals(mode, x, y, eps, masks):
    """Sz-based quantities. masks: dict name -> boolean array on the cell-centre grid.
    Returns power fractions, A_eff (Sz definition) and per-mask <Sz^2> weights for gamma."""
    dA = cell_weights(x, y)
    Ex, Ey = mode['Ex'], mode['Ey']            # already on cell centres in EMpy
    Hx, Hy = centre(mode['Hx']), centre(mode['Hy'])
    if Ex.shape != Hx.shape:
        Ex, Ey = centre(Ex), centre(Ey)
    Sz = 0.5*np.real(Ex*np.conj(Hy) - Ey*np.conj(Hx))
    P = np.sum(Sz*dA)
    out = dict(P=P, Aeff=P**2/np.sum(Sz**2*dA))
    for name, mk in masks.items():
        out['frac_' + name] = np.sum(Sz[mk]*dA[mk])/P
        out['S2_' + name] = np.sum(Sz[mk]**2*dA[mk])/P**2
    out['Sz'] = Sz
    return out
