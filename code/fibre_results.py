"""
Production runs for the liquid-core PCF (all numbers quoted in the paper come from the JSON files written here).
usage: python3 fibre_results.py <stage>
stages: design | sweep | sens | conv | modes | rings | fsm | kaggle
"""
import sys, json, time
import numpy as np
from pcf7 import *

# ---------------- design ----------------
P, D, RC, TW, NCLAD = 1.66, 1.58, 2.45, 0.08, 5      # um
LAM0 = 1.064                                          # um
DN_S = 0.016                                          # serum index offset above water
GRID = dict(h_fine=0.025, h_coarse=0.15, jacket=1.0, pml=2.5, smax=3.0)
OUT = '../out/'

def n_a(lam, dn=DN_S):
    return n_water(lam) + dn

def solve(lam, dn=DN_S, pol='x', nclad=NCLAD, grid=None, neigs=8, **kw):
    g = dict(GRID); g.update(kw)
    na = n_a(lam, dn)
    return core_mode7(lam, P, D, RC, TW, nclad, na, pol=pol, grid=grid, neigs=neigs, n_target=na - 0.008, **g)

def derivs(lam0, dn=DN_S, rel=0.02, grid=None, nclad=NCLAD):
    """beta1, beta2, beta3 at lam0 from a 5-point stencil in angular frequency"""
    w0 = 2*np.pi*c0/(lam0*1e-6); dw = rel*w0
    bs = []; mid = None
    for q in range(-2, 3):
        wq = w0 + q*dw; lam = 2*np.pi*c0/wq*1e6
        r = solve(lam, dn, grid=grid, nclad=nclad)
        bs.append(r['neff'].real*wq/c0)
        if q == 0: mid = r
    bs = np.array(bs)
    b1 = (bs[0] - 8*bs[1] + 8*bs[3] - bs[4])/(12*dw)
    b2 = (-bs[0] + 16*bs[1] - 30*bs[2] + 16*bs[3] - bs[4])/(12*dw**2)
    b3 = (bs[4] - 2*bs[3] + 2*bs[1] - bs[0])/(2*dw**3)
    return mid, b1, b2, b3

def pack(r):
    return dict(neff=float(r['neff'].real), neff_imag=float(abs(r['neff'].imag)), core=float(r['frac_core']),
                silica=float(r['frac_silica']), air=float(r['frac_air']), Aeff=float(r['Aeff']),
                S2_core=float(r['S2_core']), S2_silica=float(r['S2_silica']))

stage = sys.argv[1]
t0 = time.time()
grid = grid7(P, D, RC, NCLAD, **GRID)

if stage == 'design':
    r, b1, b2, b3 = derivs(LAM0, grid=grid)
    out = pack(r); out.update(lam0=LAM0, n_water=float(n_water(LAM0)), n_silica=float(n_silica(LAM0)), n_a=float(n_a(LAM0)),
                              ng=float(b1*c0), beta2=float(b2), beta3=float(b3), ngrid=len(grid))
    # y-polarised partner (degeneracy check)
    ry = solve(LAM0, pol='y', grid=grid)
    out['neff_y'] = float(ry['neff'].real)
    json.dump(out, open(OUT + 'fibre_design.json', 'w'), indent=1)
    # save field for plotting
    np.savez_compressed(OUT + 'fibre_mode.npz', xc=r['xc'].real, yc=r['yc'].real, Sz=r['I']['Sz'], eps=np.real(r['eps']),
                        Ex=r['mode']['Ex'], Ey=r['mode']['Ey'])
    print(json.dumps(out, indent=1))

elif stage == 'sweep':
    # uniform grid in angular frequency between 0.86 and 1.34 um
    w_lo, w_hi = 2*np.pi*c0/1.34e-6, 2*np.pi*c0/0.86e-6
    ws = np.linspace(w_lo, w_hi, 25)
    rows = []
    for wq in ws:
        lam = 2*np.pi*c0/wq*1e6
        r = solve(lam, grid=grid)
        row = pack(r); row.update(lam=float(lam), omega=float(wq), n_a=float(n_a(lam)))
        rows.append(row); print('%.4f um  neff=%.6f core=%.4f Aeff=%.3f  (%.0fs)' % (lam, row['neff'], row['core'], row['Aeff'], time.time() - t0), flush=True)
        json.dump(rows, open(OUT + 'fibre_sweep.json', 'w'), indent=1)

elif stage == 'sens':
    rows = []
    for dn in (0.0, 0.004, 0.008, 0.012, 0.0155, 0.016, 0.0165, 0.020, 0.024):
        r = solve(LAM0, dn, grid=grid)
        row = pack(r); row.update(dn=dn, n_a=float(n_a(LAM0, dn)))
        rows.append(row); print('dn=%.4f n_a=%.5f neff=%.7f core=%.4f Aeff=%.3f (%.0fs)' % (dn, row['n_a'], row['neff'], row['core'], row['Aeff'], time.time() - t0), flush=True)
    json.dump(rows, open(OUT + 'fibre_sens.json', 'w'), indent=1)
    # group index and beta2 at two neighbouring analyte indices -> d n_g / d n_a, d beta2 / d n_a
    rows2 = []
    for dn in (0.012, 0.020):
        r, b1, b2, b3 = derivs(LAM0, dn, grid=grid)
        rows2.append(dict(dn=dn, ng=float(b1*c0), beta2=float(b2), beta3=float(b3), Aeff=float(r['Aeff']), neff=float(r['neff'].real)))
        print(rows2[-1], flush=True)
    json.dump(rows2, open(OUT + 'fibre_sens_ng.json', 'w'), indent=1)

elif stage == 'conv':
    rows = []
    for hf in (0.05, 0.04, 0.03, 0.025, 0.02, 0.016):
        g = dict(GRID); g['h_fine'] = hf
        gr = grid7(P, D, RC, NCLAD, **g)
        r = solve(LAM0, grid=gr)
        row = pack(r); row.update(h=hf, ngrid=len(gr)); rows.append(row)
        print('h=%.3f N=%d neff=%.7f core=%.5f Aeff=%.4f (%.0fs)' % (hf, len(gr), row['neff'], row['core'], row['Aeff'], time.time() - t0), flush=True)
        json.dump(rows, open(OUT + 'fibre_conv.json', 'w'), indent=1)

elif stage == 'modes':
    # first few core modes in each symmetry class (to document that the liquid core is multimode)
    from fibre_tools import solve_modes, mode_integrals
    na = n_a(LAM0); ns = n_silica(LAM0)
    res = {}
    for bc in ('0A0A', '0S0S', '0A0S', '0S0A'):
        modes, eps, (xc, yc) = solve_modes(LAM0, grid, grid.copy(), shapes7(P, D, RC, TW, NCLAD, na, ns), ns**2, bc, neigs=40, n_target=na - 0.03, ss=6)
        XC, YC = np.meshgrid(xc.real, yc.real, indexing='ij'); core = XC**2 + YC**2 <= RC**2
        lst = []
        for mo in modes:
            I = mode_integrals(mo, grid, grid, eps, dict(core=core))
            if I['frac_core'] > 0.5:
                lst.append((float(mo['neff'].real), float(abs(mo['neff'].imag)), float(I['frac_core'])))
        res[bc] = sorted(lst, reverse=True)
        print(bc, res[bc], '(%.0fs)' % (time.time() - t0), flush=True)
    json.dump(res, open(OUT + 'fibre_modes.json', 'w'), indent=1)

elif stage == 'rings':
    rows = []
    for nc in (1, 2, 3):
        g = dict(GRID)
        gr = grid7(P, D, RC, nc, **g)
        r = solve(LAM0, nclad=nc, grid=gr, neigs=12)
        row = pack(r); row.update(nclad=nc); rows.append(row)
        print('nclad=%d neff=%.7f Im=%.3e core=%.5f (%.0fs)' % (nc, row['neff'], row['neff_imag'], row['core'], time.time() - t0), flush=True)
    json.dump(rows, open(OUT + 'fibre_rings.json', 'w'), indent=1)

elif stage == 'fsm':
    from fsm import n_fsm
    out = {}
    for lam in (0.90, 1.064, 1.30):
        v, both = n_fsm(lam, P, D, h=0.008)
        out[str(lam)] = float(v); print(lam, v, both, flush=True)
    json.dump(out, open(OUT + 'fibre_fsm.json', 'w'), indent=1)
elif stage == 'kaggle':
    # mode solver evaluated at the 18 distinct blood-glucose levels of the Kaggle "Diabetes prediction dataset"
    # (out/kaggle_levels.json is written by kaggle_summary.py); reference level C_REF has the nominal serum index.
    DNDC, C_REF = 1.4e-6, 100.0                        # RIU per mg/dL (Martini et al. 2009); mg/dL
    K = json.load(open(OUT + 'kaggle_levels.json'))
    rows = []
    for C in K['levels']:
        dn = DN_S + DNDC*(C - C_REF)
        r = solve(LAM0, dn, grid=grid)
        row = pack(r); row.update(C_mgdl=C, dn=dn, n_a=float(n_a(LAM0, dn)))
        rows.append(row); print('C=%3d mg/dL n_a=%.7f neff=%.9f core=%.5f Aeff=%.4f (%.0fs)' % (C, row['n_a'], row['neff'], row['core'], row['Aeff'], time.time() - t0), flush=True)
        json.dump(rows, open(OUT + 'fibre_kaggle.json', 'w'), indent=1)
print('stage %s done in %.0f s' % (stage, time.time() - t0))
