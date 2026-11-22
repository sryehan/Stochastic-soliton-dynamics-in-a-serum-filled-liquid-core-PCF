"""Collect every number quoted in the paper from the raw solver output."""
import json, numpy as np
c0 = 299792458.0
O = '../out/'
d = json.load(open(O + 'fibre_design.json')); sw = json.load(open(O + 'fibre_sweep.json')); se = json.load(open(O + 'fibre_sens.json'))
sg = json.load(open(O + 'fibre_sens_ng.json')); cv = json.load(open(O + 'fibre_conv.json')); rg = json.load(open(O + 'fibre_rings.json'))
fs = json.load(open(O + 'fibre_fsm.json')); md = json.load(open(O + 'fibre_modes.json'))
N = {}
lam0 = d['lam0']*1e-6; k0 = 2*np.pi/lam0
# material constants (literature)
n2_w, n2_s = 4.1e-20, 2.6e-20          # m^2/W  (water: Boyd; silica: Agrawal)
alpha_w = 14.2                         # 1/m at 1064 nm (Peterman et al. 2003)
dndC = 1.4e-6                          # RIU per mg/dL (Martini et al. 2009)
dndT = -1.0e-4                         # 1/K
N.update(lam0_um=d['lam0'], n_water=d['n_water'], n_silica=d['n_silica'], n_a=d['n_a'], neff=d['neff'], neff_y=d['neff_y'],
         birefringence_numerical=abs(d['neff'] - d['neff_y']), core_fraction=d['core'], silica_fraction=d['silica'], air_fraction=d['air'],
         Aeff_um2=d['Aeff'], ng=d['ng'], beta2_ps2_km=d['beta2']*1e27, beta3_ps3_km=d['beta3']*1e39, n_fsm=fs['1.064'])
N['D_ps_nm_km'] = -2*np.pi*c0/lam0**2*d['beta2']*1e6
gamma = k0*(n2_w*d['S2_core'] + n2_s*d['S2_silica'])*1e12
N['gamma_W_m'] = gamma; N['gamma_simple'] = k0*n2_w/(d['Aeff']*1e-12); N['gamma_silica_share'] = n2_s*d['S2_silica']/(n2_w*d['S2_core'] + n2_s*d['S2_silica'])
# sensitivity
dn = np.array([r['dn'] for r in se]); ne = np.array([r['neff'] for r in se]); Ae = np.array([r['Aeff'] for r in se])
p = np.polyfit(dn, ne, 1); N['Sn_linear_fit'] = p[0]; N['Sn_fit_maxres'] = float(np.max(np.abs(np.polyval(p, dn) - ne)))
i = list(dn).index(0.016); N['Sn_local'] = (ne[i+1] - ne[i-1])/(dn[i+1] - dn[i-1])
N['dlnAeff_dna'] = (Ae[i+1] - Ae[i-1])/(dn[i+1] - dn[i-1])/Ae[i]
N['dng_dna'] = (sg[1]['ng'] - sg[0]['ng'])/(sg[1]['dn'] - sg[0]['dn'])
N['dlnbeta2_dna'] = (sg[1]['beta2'] - sg[0]['beta2'])/(sg[1]['dn'] - sg[0]['dn'])/d['beta2']
N['dlnP0_dna'] = N['dlnbeta2_dna'] + N['dlnAeff_dna']
alpha = N['Sn_local']*alpha_w; N['alpha_mode_1_m'] = alpha; N['alpha_dB_cm'] = 4.343*alpha/100; N['L_abs_cm'] = 100/alpha
# confinement loss
N['rings'] = [dict(nclad=r['nclad'], neff_imag=r['neff_imag'], loss_dB_m=8.686*k0*r['neff_imag']) for r in rg]
N['conv'] = [dict(h=r['h'], neff=r['neff'], Aeff=r['Aeff']) for r in cv]
N['conv_neff_err_h025'] = abs(cv[3]['neff'] - cv[-1]['neff']); N['conv_Aeff_relerr_h025'] = abs(cv[3]['Aeff'] - cv[-1]['Aeff'])/cv[-1]['Aeff']
# dispersion curve
om = np.array([r['omega'] for r in sw]); lam = np.array([r['lam'] for r in sw]); be = np.array([r['neff'] for r in sw])*om/c0
h = om[1] - om[0]
b2 = np.full_like(be, np.nan); b2[2:-2] = (-be[:-4] + 16*be[1:-3] - 30*be[2:-2] + 16*be[3:-1] - be[4:])/(12*h*h)
b1 = np.full_like(be, np.nan); b1[2:-2] = (be[:-4] - 8*be[1:-3] + 8*be[3:-1] - be[4:])/(12*h)
ok = ~np.isnan(b2)
# zero-dispersion wavelength by linear interpolation of beta2(omega)
zz = None
for j in range(len(om) - 1):
    if ok[j] and ok[j+1] and b2[j]*b2[j+1] < 0:
        wz = om[j] - b2[j]*(om[j+1] - om[j])/(b2[j+1] - b2[j]); zz = 2*np.pi*c0/wz*1e6
N['ZDW_um'] = zz
N['sweep'] = [dict(lam=float(l), neff=float(r['neff']), core=r['core'], Aeff=r['Aeff'], beta2_ps2_km=(None if np.isnan(x) else float(x*1e27)), ng=(None if np.isnan(y) else float(y*c0)))
              for l, r, x, y in zip(lam, sw, b2, b1)]
N['beta2_interp_at_lam0'] = float(np.interp(2*np.pi*c0/lam0, om[ok], b2[ok])*1e27)
# solitons
tab = []
for T0 in (20e-15, 30e-15, 50e-15):
    LD = T0**2/abs(d['beta2']); P0 = 1/(gamma*LD)
    tab.append(dict(T0_fs=T0*1e15, FWHM_fs=1.7627*T0*1e15, LD_cm=LD*100, P0_W=P0, E_pJ=2*P0*T0*1e12, Gamma=alpha*LD/2, alphaLD=alpha*LD,
                    delta3=d['beta3']/(6*abs(d['beta2'])*T0), s=lam0/(2*np.pi*c0*T0), I_GW_cm2=P0/(d['Aeff']*1e-12)*1e-13,
                    LD3_cm=T0**3/d['beta3']*100, z0_cm=np.pi/2*LD*100, FOM=1/(alpha*LD)))
N['solitons'] = tab
# sensing, nominal: T0 = 30 fs, L = 5 cm
T0 = 30e-15; LD = T0**2/abs(d['beta2']); L = 0.05
N.update(L_cm=L*100, L_over_LD=L/LD, transmission=float(np.exp(-alpha*L)), S_phi_rad_per_RIU=k0*L*N['Sn_local'],
         dphi_per_mgdl=k0*L*N['Sn_local']*dndC, S_tau_ps_per_RIU=L/c0*N['dng_dna']*1e12, dtau_fs_per_mgdl=L/c0*N['dng_dna']*dndC*1e15,
         kappa_per_mgdl=k0*N['Sn_local']*dndC*LD, phi_NL=L/(2*LD), dlnP0_per_mgdl=N['dlnP0_dna']*dndC,
         glucose_equiv_per_K=abs(dndT)/dndC, dphi_per_K=k0*L*N['Sn_local']*abs(dndT), dn_300mgdl=300*dndC)
N['LoD_mgdl'] = {str(v): 3*v*np.sqrt(L/LD)/(k0*L*N['Sn_local']*dndC) for v in (0.01, 0.05, 0.1, 0.2, 0.3)}
N['LoD_RIU'] = {str(v): 3*v*np.sqrt(L/LD)/(k0*L*N['Sn_local']) for v in (0.01, 0.05, 0.1, 0.2, 0.3)}
# thermal-noise illustration: delta T rms, correlation length lc
for dT, lc in ((1e-3, 1e-3), (1e-2, 1e-3)):
    sb = k0*N['Sn_local']*abs(dndT)*dT; Dphi = 2*sb*sb*lc
    N['varrho_dT%g' % dT] = float(np.sqrt(Dphi*LD))
N['modes'] = {k: [(round(a, 5), round(c_, 3)) for a, b_, c_ in v[:4]] for k, v in md.items()}
N['consts'] = dict(n2_w=n2_w, n2_s=n2_s, alpha_w=alpha_w, dndC=dndC, dndT=dndT)
json.dump(N, open(O + 'numbers.json', 'w'), indent=1, default=float)
for k, v in N.items():
    if k not in ('sweep',): print(k, ':', v)
