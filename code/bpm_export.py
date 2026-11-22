"""Converts the output of matlab/m6_beam_propagation.m (bpm_results.mat) into the compact files used by fig_bpm.py:
out/bpm_results.json (numbers and curves) and out/bpm_snapshots.npz (index maps and intensity pictures).

Usage:  python bpm_export.py [path to bpm_results.mat, default ../matlab/bpm_results.mat]
"""
import sys, json
import numpy as np, scipy.io as sio
from fsm import n_fsm

path = sys.argv[1] if len(sys.argv) > 1 else '../matlab/bpm_results.mat'
m = sio.loadmat(path)
x = m['x'].ravel(); sel = np.abs(x) <= 11.5
sc = lambda k: float(np.asarray(m[k]).ravel()[0])
res = [m['res'][0, q] for q in range(3)]
get = lambda s, k: np.asarray(s[k][0, 0])
out = dict(designs=['A', 'B', 'C'], Z_um=sc('Z'), Zg_um=sc('Zg'), dz_um=sc('dz'), w0_um=sc('w0'), pixel_um=float(x[1] - x[0]), N=int(len(x)),
           z=get(res[0], 'z').ravel().round(4).tolist(),
           Pcore=[get(r, 'Pcore').ravel().tolist() for r in res], Ptot=[get(r, 'Ptot').ravel().tolist() for r in res],
           snap_z=get(res[0], 'snap_z').ravel().tolist(),
           decay_A_per_um=sc('aA'), decay_A_exact_per_um=sc('aEx'), decay_A_MS_scalar=sc('aMS_s'), decay_A_MS_EH11=sc('aMS_v'),
           decay_B_window_per_um=sc('aB'),
           mode_C=dict(neff=sc('neffC'), core_fraction=sc('coreC'), Aeff_um2=sc('AeffC'), launch_overlap=sc('eta')),
           glucose=dict(levels_mgdl=np.asarray(m['Cg']).ravel().tolist(), z=m['zg'].ravel().round(4).tolist(), dphi=m['dphi'].ravel().tolist(),
                        slope_rad_per_um=float(m['pg'].ravel()[0]), Sn=sc('Sn'),
                        Ptot_end=[float(get(m['og'][0, q], 'Ptot').ravel()[-1]) for q in range(2)]))
zz = np.array(out['z'])
for q, d in enumerate('ABC'):
    out['Pcore_' + d] = {str(int(v)): float(np.interp(v, zz, out['Pcore'][q])) for v in (10, 20, 50, 100, 150, 200) if v <= out['Z_um']}
    out['Ptot_' + d] = {str(int(v)): float(np.interp(v, zz, out['Ptot'][q])) for v in (10, 20, 50, 100, 150, 200) if v <= out['Z_um']}
half = zz >= out['Z_um']/2
out['nfsm_B'] = float(n_fsm(1.064, 4.46, 2.23, h=0.02)[0])          # cladding index of design B (fundamental space-filling mode)
out['C_core_mean_last_half'] = float(np.mean(np.array(out['Pcore'][2])[half])); out['C_core_std_last_half'] = float(np.std(np.array(out['Pcore'][2])[half]))
json.dump(out, open('../out/bpm_results.json', 'w'), indent=1)
maps = np.stack([np.sqrt(np.asarray(m['maps'][0, q]))[np.ix_(sel, sel)] for q in range(3)]).astype(np.float32)
snaps = np.stack([[10*np.log10(np.asarray(get(r, 'snap')[0, j])[np.ix_(sel, sel)] + 1e-12) for j in range(3)] for r in res])
snaps_u8 = np.clip(np.round((snaps + 51.0)*5), 0, 255).astype(np.uint8)        # dB = value/5 - 51, 0.2 dB steps
np.savez_compressed('../out/bpm_snapshots.npz', x=x[sel].astype(np.float32), index=maps, intensity_dB_u8=snaps_u8)
print('wrote ../out/bpm_results.json and ../out/bpm_snapshots.npz')
print('A: decay %.4f per um (exact leaky mode %.4f); B: core %.3f at 50 um, %.3f at 200 um; C: core %.4f +- %.4f, window %.4f at 200 um'
      % (out['decay_A_per_um'], out['decay_A_exact_per_um'], out['Pcore_B']['50'], out['Pcore_B']['200'],
         out['C_core_mean_last_half'], out['C_core_std_last_half'], out['Ptot_C']['200']))
