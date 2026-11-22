"""Collect the data of Figs. 5-12 of the paper into one MATLAB file, out/figure_data.mat, for matlab/paper_figures.m.
Everything comes from the files in out/ written by the computations; nothing is recomputed here."""
import json, numpy as np, scipy.io as sio
from fibre_tools import n_water

O = '../out/'
J = lambda f: json.load(open(O + f))
D = {}

# Fig. 5: stochastic verification
P, E, W = J('stoch_path.json'), J('stoch_ens.json'), J('pathwise.json')
D['f5'] = dict(t=np.array(P['t']), re_det=np.array(P['re_det']), re_num=np.array(P['re_num']), re_ex=np.array(P['re_ex']),
               I_num=np.array(P['I_num']), bright=np.array(W['bright']), dark=np.array(W['dark']),
               bright_order=W['bright_order'], dark_order=W['dark_order'], te=np.array(E['t']),
               rho=np.array([0.1, 0.2, 0.3]),
               meanfield=np.array([E[v]['meanfield'] for v in ('0.1', '0.2', '0.3')]),
               phvar=np.array([E[v]['phvar'] for v in ('0.1', '0.2', '0.3')]))

# Fig. 6: fundamental mode (one quadrant of S_z)
Z = np.load(O + 'fibre_mode.npz')
D['f6'] = dict(xc=Z['xc'], yc=Z['yc'], Sz=Z['Sz']/Z['Sz'].max())

# Fig. 7: dispersion and effective area
N = J('numbers.json'); sw = N['sweep']
lam = np.array([r['lam'] for r in sw])
nan = lambda v: np.nan if v is None else v
D['f7'] = dict(lam=lam, neff=np.array([r['neff'] for r in sw]), ng=np.array([nan(r['ng']) for r in sw]),
               beta2=np.array([nan(r['beta2_ps2_km']) for r in sw]), Aeff=np.array([r['Aeff'] for r in sw]),
               n_a=n_water(lam) + 0.016, ZDW=N['ZDW_um'], beta2_design=N['beta2_ps2_km'], Aeff_design=N['Aeff_um2'])

# Fig. 8: propagation with loss and third-order dispersion
R = J('fibre_prop.json'); f8 = {}
for key in ('30', '20'):
    r = R[key]; f0 = r['ideal']['fwhm'][0]
    f8['T%s' % key] = dict(z_cm=np.array(r['ideal']['xi'])*r['LD_cm'], linear=np.array(r['linear']['fwhm'])/f0,
                           soliton=np.array(r['loss_tod']['fwhm'])/f0, ideal=np.array(r['ideal']['fwhm'])/f0)
r = R['30']
f8['snap'] = dict(T_fs=np.array(r['loss_tod']['x_snap'])*30.0, soliton=np.array(r['loss_tod']['I_snap']), linear=np.array(r['linear']['I_snap']))
D['f8'] = f8

# Fig. 9: beam propagation (snapshots stored as dB, 0 ... -51 dB, in steps of 0.2 dB)
B = J('bpm_results.json'); S = np.load(O + 'bpm_snapshots.npz')
D['f9'] = dict(x=S['x'], index=S['index'].astype(np.float32), intensity_dB=(S['intensity_dB_u8'].astype(np.float32)/5 - 51.0),
               snap_z=np.array(B['snap_z']), z=np.array(B['z']), Pcore=np.array(B['Pcore']), Ptot=np.array(B['Ptot']),
               decay_A_exact=B['decay_A_exact_per_um'], ng=N['ng'],
               gz=np.array(B['glucose']['z']), gdphi=np.array(B['glucose']['dphi']), gslope=B['glucose']['slope_rad_per_um'])

# Fig. 10: side view of the beam propagation (bpm_movie.py)
V = np.load(O + 'bpm_side_view.npz')
D['f10'] = dict(side_dB=(10*np.log10(V['side'] + 1e-12)).astype(np.float32), xs=V['xs'], zs=V['zs'], Rc=float(V['Rc']))

# Fig. 11: sensing
se = J('fibre_sens.json')
D['f11'] = dict(n_a=np.array([r['n_a'] for r in se]), neff=np.array([r['neff'] for r in se]),
                dphi_per_mgdl=N['dphi_per_mgdl'], dtau_fs_per_mgdl=N['dtau_fs_per_mgdl'],
                LD_m=N['solitons'][1]['LD_cm']/100, Sn=N['Sn_local'], dndC=N['consts']['dndC'], lam0_m=1.064e-6)

# Fig. 12: data set
A = J('kaggle_analysis.json')
D['f12'] = dict(levels=np.array(A['levels']), counts=np.array(A['counts']), counts_diabetes=np.array(A['counts_diabetes']),
                neff=np.array(A['neff']), sigT=np.array(A['curve_sigT']),
                within_rho01=np.array(A['curves']['0.1']), within_rho10=np.array(A['curves']['1.0']))

sio.savemat(O + 'figure_data.mat', D, do_compression=True, long_field_names=True)
print('wrote out/figure_data.mat')
