"""Videos of the beam travelling along the three liquid-core fibers (Python version of matlab/m7_propagation_movie.m).

The same Gaussian beam (1064 nm) is launched into A (liquid core in solid silica), B (conventional photonic-crystal
cladding, d/Lambda = 0.5) and C (this work, d/Lambda = 0.952) and propagated over 200 um with the scalar
split-step beam propagation method of matlab/bpm_run.m (same grid, step, absorber and launched beam).

Writes to ../videos/ (H.264 .mp4, 30 frames/s; needs ffmpeg on the PATH):
  A_solid_silica.mp4, B_conventional_cladding.mp4, C_this_work.mp4   1080 x 1080
  ABC_comparison.mp4                                                 1920 x 720
  side_view.png                                                      intensity on y = 0 against distance
Run time: about 10 min.   Usage:  python bpm_movie.py  [--fast]
"""
import os, sys, subprocess, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
try:
    from scipy import fft as sfft
    fft2 = lambda a: sfft.fft2(a, workers=-1); ifft2 = lambda a: sfft.ifft2(a, workers=-1)
except ImportError:
    fft2, ifft2 = np.fft.fft2, np.fft.ifft2

FAST = '--fast' in sys.argv
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'videos')
os.makedirs(OUT, exist_ok=True)

# ---------------- parameters (as in matlab/m6_beam_propagation.m) ----------------
lam0 = 1.064; k0 = 2*np.pi/lam0
n_si, n_w = 1.449631, 1.323917; n_a = n_w + 0.016; n0 = 1.332
N, dx = 640, 0.05; x = (np.arange(N) - N//2)*dx
dz = 0.05 if FAST else 0.025
Z, w0, Rc, rabs = 200.0, 1.9, 2.45, 11.0
c0, n_g = 299.792458, 1.3633
dzf, fps, hold_end = 1.0, 30, 2.0
TITLES = ['A: liquid core in solid silica', 'B: conventional cladding, $d/\\Lambda=0.5$', 'C: this work, $d/\\Lambda=0.952$']
NAMES = ['A_solid_silica', 'B_conventional_cladding', 'C_this_work']
FOOTER = 'scalar beam propagation, $\\lambda$ = 1064 nm, serum-filled core; colour: intensity in dB below the launched peak'
CM = LinearSegmentedColormap.from_list('bpm', np.array([[255, 255, 255], [207, 225, 247], [127, 176, 234],
                                                         [42, 120, 214], [18, 63, 122], [7, 29, 58]])/255)

def lattice(P, r1, r2):
    xy = []
    for i in range(-r2, r2 + 1):
        for j in range(-r2, r2 + 1):
            r = max(abs(i), abs(j), abs(i + j))
            if r1 <= r <= r2: xy.append(((i + j/2)*P, j*np.sqrt(3)/2*P))
    return xy

def index_map(design, s=4):
    """n^2(x, y) of matlab/pcf_index_map.m, averaged over s x s sub-pixels."""
    xf = x[0] - dx/2 + (np.arange(N*s) + 0.5)*dx/s
    XF, YF = np.meshgrid(xf, xf)
    ef = np.full((N*s, N*s), n_si**2)
    if design == 'A': holes, rh, rwall = [], 0, 0
    elif design == 'B': holes, rh, rwall = lattice(4.46, 1, 2), 2.23/2, 0
    else: holes, rh, rwall = lattice(1.66, 2, 6), 1.58/2, Rc + 0.08
    for (hx, hy) in holes: ef[(XF - hx)**2 + (YF - hy)**2 <= rh**2] = 1.0
    if rwall: ef[XF**2 + YF**2 <= rwall**2] = n_si**2
    ef[XF**2 + YF**2 <= Rc**2] = n_a**2
    return ef.reshape(N, s, N, s).mean(axis=(1, 3))

X, Y = np.meshgrid(x, x); R = np.hypot(X, Y)
kx = 2*np.pi*np.fft.fftfreq(N, dx); KX, KY = np.meshgrid(kx, kx); K2 = KX**2 + KY**2
D = np.exp(-1j*K2/(2*k0*n0)*dz); Dh = np.exp(-1j*K2/(2*k0*n0)*dz/2)
sig = 1.0*np.maximum(0, (R - rabs)/(x.max() - rabs))**2
psi0 = np.exp(-R**2/w0**2); I0 = np.max(np.abs(psi0)**2)
sel = np.where(np.abs(x) <= 11.5)[0]; xs = x[sel]; j0 = np.argmin(np.abs(x))
ext = [xs[0] - dx/2, xs[-1] + dx/2, xs[0] - dx/2, xs[-1] + dx/2]
maps = [index_map(d) for d in 'ABC']
V = [np.exp(1j*k0*(m - n0**2)/(2*n0)*dz - sig*dz) for m in maps]
psi = [ifft2(Dh*fft2(psi0)) for _ in range(3)]

# ---------------- figures ----------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13})
dB = lambda I: 10*np.log10(I + 1e-12)
def panel(ax, q, I):
    im = ax.imshow(dB(I), extent=ext, origin='lower', cmap=CM, vmin=-30, vmax=0, interpolation='bilinear')
    ax.contour(xs, xs, np.sqrt(maps[q][np.ix_(sel, sel)]), [1.2, 1.39], colors=[(0.6, 0.6, 0.6)], linewidths=0.6)
    ax.set_xlabel('x (µm)'); ax.set_ylabel('y (µm)')
    return im

def ffmpeg(name, w, h):
    return subprocess.Popen(['ffmpeg', '-y', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgba', '-s', '%dx%d' % (w, h),
                             '-r', str(fps), '-i', '-', '-c:v', 'libx264', '-preset', 'slow', '-crf', '16', '-pix_fmt', 'yuv420p',
                             '-movflags', '+faststart', os.path.join(OUT, name + '.mp4')], stdin=subprocess.PIPE)

I_start = np.abs(psi0[np.ix_(sel, sel)])**2/I0
figW = plt.figure(figsize=(16, 6), dpi=120)
axW = [figW.add_axes([0.04 + q*0.295, 0.15, 0.26, 0.66]) for q in range(3)]
imW = [panel(axW[q], q, I_start) for q in range(3)]
for q in range(3):
    axW[q].set_title(TITLES[q], fontsize=15)
    if q: axW[q].set_ylabel('')
cb = figW.colorbar(imW[2], cax=figW.add_axes([0.925, 0.15, 0.011, 0.66])); cb.set_label('dB')
tW = figW.text(0.5, 0.925, '', ha='center', fontsize=19, fontweight='bold')
figW.text(0.5, 0.025, FOOTER, ha='center', fontsize=11, color='0.35')

figS = plt.figure(figsize=(9, 9), dpi=120)
axS = figS.add_axes([0.11, 0.12, 0.72, 0.72])
imS = axS.imshow(dB(I_start), extent=ext, origin='lower', cmap=CM, vmin=-30, vmax=0, interpolation='bilinear')
axS.set_xlabel('x (µm)'); axS.set_ylabel('y (µm)')
contS = []
cbS = figS.colorbar(imS, cax=figS.add_axes([0.855, 0.12, 0.022, 0.72])); cbS.set_label('intensity (dB)')
hT = axS.set_title('', fontsize=17)
tS = figS.text(0.5, 0.925, '', ha='center', fontsize=19, fontweight='bold')
figS.text(0.5, 0.025, FOOTER.replace('; ', ';\n'), ha='center', fontsize=10, color='0.35')

def grab(fig):
    fig.canvas.draw(); b = np.asarray(fig.canvas.buffer_rgba()); return b.tobytes(), b.shape

_, shW = grab(figW); _, shS = grab(figS)
pW = ffmpeg('ABC_comparison', shW[1], shW[0]); pS = [ffmpeg(n, shS[1], shS[0]) for n in NAMES]

# ---------------- propagate A, B, C together ----------------
every = int(round(dzf/dz)); nst = int(round(Z/dz)); everyS = int(round(0.5/dz))
zs = np.arange(0, Z + 1e-9, 0.5); side = np.zeros((len(sel), len(zs), 3)); ks = 0
Icur = [I_start.copy() for _ in range(3)]
for q in range(3): side[:, 0, q] = np.abs(psi0[j0, sel])**2/I0
t0 = time.time()
for it in range(nst + 1):
    if it:
        rec, recS = it % every == 0, it % everyS == 0
        for q in range(3):
            psi[q] *= V[q]
            if rec or recS:
                I = np.abs(ifft2(Dh*fft2(psi[q])))**2/I0
                if recS: side[:, ks + 1, q] = I[j0, sel]
                if rec: Icur[q] = I[np.ix_(sel, sel)]
            psi[q] = ifft2(D*fft2(psi[q]))
        if recS: ks += 1
        if not rec: continue
    z = it*dz; lbl = 'z = %5.1f µm      travel time  t = %.3f ps' % (z, z*n_g/c0)
    nrep = int(round(hold_end*fps)) if it == nst else 1
    for q in range(3): imW[q].set_data(dB(Icur[q]))
    tW.set_text(lbl); buf, _ = grab(figW)
    for _ in range(nrep): pW.stdin.write(buf)
    for q in range(3):
        imS.set_data(dB(Icur[q]))
        for c in contS: c.remove()
        cs = axS.contour(xs, xs, np.sqrt(maps[q][np.ix_(sel, sel)]), [1.2, 1.39], colors=[(0.6, 0.6, 0.6)], linewidths=0.6)
        contS = [cs]; hT.set_text(TITLES[q]); tS.set_text(lbl)
        buf, _ = grab(figS)
        for _ in range(nrep): pS[q].stdin.write(buf)
    if abs(z % 20) < 1e-9: print('z = %3.0f um  (%.0f s)' % (z, time.time() - t0), flush=True)
for p in [pW] + pS:
    p.stdin.close(); p.wait()

# ---------------- side view ----------------
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), dpi=150, constrained_layout=True)
for q in range(3):
    im = axs[q].imshow(dB(side[:, :, q]), extent=[0, Z, ext[2], ext[3]], origin='lower', aspect='auto', cmap=CM, vmin=-30, vmax=0)
    axs[q].axhline(Rc, color='k', ls=':', lw=0.8); axs[q].axhline(-Rc, color='k', ls=':', lw=0.8)
    axs[q].set_title(TITLES[q]); axs[q].set_xlabel('distance z (µm)')
axs[0].set_ylabel('x (µm), on the line y = 0')
fig.colorbar(im, ax=axs, label='intensity (dB)', shrink=0.9)
fig.savefig(os.path.join(OUT, 'side_view.png'))
np.savez_compressed(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'out', 'bpm_side_view.npz'),
                    side=side.astype(np.float32), xs=xs, zs=zs, Rc=Rc)
print('done: videos and side view in', os.path.abspath(OUT), '(%.0f s)' % (time.time() - t0))
