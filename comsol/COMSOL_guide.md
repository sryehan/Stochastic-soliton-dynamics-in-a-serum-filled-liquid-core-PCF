# Building the liquid-core PCF in COMSOL Multiphysics

This folder lets a co-author repeat, in COMSOL, the mode calculation of Section 6 of the paper. The results in the
paper come from a full-vector finite-difference solver that is validated against an exact analytical solution and a
published multipole benchmark (Table 5). COMSOL is therefore a **cross-check**, not the source of any number.

| File | What it is |
|---|---|
| `lcpcf_model.m` | Complete model for LiveLink for MATLAB. **Never executed by its author** — no COMSOL licence was available. |
| `lcpcf_cross_section.dxf` | The cross-section as circles, for `Import` in a 2D geometry (base COMSOL reads 2D DXF). |
| `lcpcf_holes.csv` | All 120 hole centres, with their ring index and distance from the axis. |
| `lcpcf_sector.csv` | The 20 hole centres of one 60° sector, for the quickest manual build. |

Two routes. With LiveLink, run `lcpcf_model.m`; if a line fails, fix that one step using Section 2 below. Without
LiveLink, build the model by hand from Section 2. Either way, compare against the checklist in Section 4.

---

## 1. Parameters

Put these in **Global Definitions > Parameters**.

| Name | Expression | Meaning |
|---|---|---|
| `lam0` | `1.064[um]` | wavelength |
| `f0` | `c_const/lam0` | mode-analysis frequency |
| `Lambda` | `1.66[um]` | pitch of the triangular lattice |
| `dh` | `1.58[um]` | air-hole diameter (`dh/Lambda` = 0.952) |
| `Rc` | `2.45[um]` | radius of the liquid core |
| `tw` | `80[nm]` | thickness of the silica core wall |
| `Rjacket` | `5*Lambda + dh/2 + 1[um]` | = 10.09 µm, outer radius of the silica region |
| `dpml` | `2.5[um]` | PML thickness |
| `n_sio2` | `1.449631` | fused silica at 1.064 µm (Malitson) |
| `n_w` | `1.323917` | water at 1.064 µm (Kedenburg) |
| `n_a` | `1.339917` | serum = `n_w + 0.016` |

For a wavelength sweep, replace the last three by analytic functions of `lam0` (Section 5.2).

## 2. Geometry (2D, length unit µm)

The cladding is a triangular lattice of air holes; the lattice sites of ring 0 and ring 1 are replaced by the
liquid core, so the holes occupy rings 2 to 6 — **120 holes**.

**2.1 Jacket and PML.** `Circle`, radius `Rjacket + dpml`, centred at the origin. Open **Layers**, add one layer
named `pml` of thickness `dpml`. This produces the silica disc with a concentric PML ring around it.

**2.2 The 120 holes.** Three ways, fastest first.

*By rotation (recommended for hand building).* Create the 20 circles of radius `dh/2` listed in
`lcpcf_sector.csv` — one 60° sector of the lattice. Then add one `Rotate` feature: select all 20 circles,
rotation angles `60 120 180 240 300`, centre the origin, and tick **Keep input objects**. That gives all 120.

*By DXF.* `Import` > `lcpcf_cross_section.dxf`. The holes are on the layer `holes`; `core`, `wall`, `jacket` and
`pml` carry the four concentric circles. Set the import length unit to µm. Check afterwards that the outermost
hole centre is at 9.96 µm from the axis.

*By LiveLink.* The loop in `lcpcf_model.m` writes the sites out of the lattice indices directly.

**2.3 Core wall.** `Circle`, radius `Rc + tw` (= 2.53 µm). Then `Difference`: **objects to add** = the 120 holes,
**objects to subtract** = this circle, with **Keep objects to subtract** ticked. This is what gives the core its
80 nm silica wall. Note that `Rc + tw` equals `2*Lambda - dh/2` exactly, so the six holes centred at 2Λ = 3.32 µm
only touch the wall, while the six centred at √3Λ = 2.875 µm are cut by it. Six tangent, six truncated — if your
geometry shows anything else, a coordinate is wrong.

**2.4 Core.** `Circle`, radius `Rc`. Build the geometry (`Build All`).

Expect **123 domains**: the core, the wall, 120 holes and the silica, plus the PML ring.

## 3. Materials, PML, physics, mesh, study

**3.1 Materials.** Add three blank materials, each with the property group **Refractive index**, and assign them
**in this order** so that later ones override earlier ones:

| Material | `n` | `k` | Selection |
|---|---|---|---|
| silica | `n_sio2` | 0 | All domains |
| air | `1` | 0 | the 120 hole domains |
| serum | `n_a` | 0 | the core domain |

The PML ring keeps the silica index; that is correct.

**3.2 PML.** `Definitions > Perfectly Matched Layer`, geometry type **Cylindrical**, assigned to the PML ring only.

**3.3 Physics.** `Electromagnetic Waves, Frequency Domain` (`ewfd`) on all domains. In **Wave Equation,
Electric**, set the electric displacement field model to **Refractive index**. No boundary condition needs to be
changed: the default perfect electric conductor sits outside the PML, which absorbs everything that reaches it.

**3.4 Mesh.** Free triangular, with these maximum element sizes:

| Region | Max element size |
|---|---|
| core wall (the 80 nm annulus) | **30 nm** |
| core and the holes of rings 2–3 | 0.15 µm |
| everything else | 0.5 µm |

The wall is the part that decides the answer: it is 80 nm thick and carries the index step that confines the
mode. If it has fewer than two or three elements across, `n_eff` will be wrong in the fourth decimal. A mapped
mesh in the PML (about 8 elements radially) is better than triangles but not essential.

**3.5 Study.** `Mode Analysis`. Mode analysis frequency `f0`. **Desired number of modes** 10, **Search for modes
around** `1.3320`. Compute.

## 4. Results, and what counts as agreement

Read the effective index as `ewfd.neff` for each mode. Add two surface integrations of the axial Poynting
component `ewfd.Poavz`:

- over **all domains**: call the values `S1 = ∫ Poavz dA` and `S2 = ∫ Poavz^2 dA`;
- over the **core only**: `Sc = ∫ Poavz dA`.

Then

- effective area `A_eff = S1^2 / S2`,
- power fraction in the serum `= Sc / S1`,
- confinement loss in dB/m `= 8.686 * (2*pi/lam0) * imag(ewfd.neff)`.

The fundamental mode is the one with the largest power fraction in the core, not simply the highest index.

### Checklist (finite-difference solver of the paper, 1.064 µm, five rings)

| Quantity | Paper | Agreement if |
|---|---|---|
| `n_eff`, fundamental mode | 1.332018 | within 1 × 10⁻⁴ |
| second polarization | 1.332014 | the two are degenerate to about 10⁻⁵ |
| power fraction in the serum | 99.19 % | within 0.2 points |
| power fraction in silica / air | 0.42 % / 0.25 % | within 0.2 points |
| `A_eff` (from `Poavz`) | 11.78 µm² | within 2 % |
| group index `n_g` | 1.3633 | within 0.001 |
| next mode group | 1.3192 – 1.3205 | 4 modes in this band |
| third group | 1.3035 – 1.3040 | present |
| confinement loss, 5 rings | below the numerical floor | ≤ 10⁻⁶ dB/m |

A COMSOL run that lands inside these bands confirms the design. The two methods use different formulations and
different meshes, so agreement in the fifth decimal of `n_eff` is not expected and not needed.

## 5. The four extra runs

**5.1 Index sensitivity `S_n`.** Run again with `n_a = 1.343917` (i.e. + 0.004) and divide the change of `n_eff`
by 0.004. Paper: **0.996**. This is the number that carries the whole sensing section, so it is the most valuable
one to confirm.

**5.2 Dispersion.** Make `n_sio2`, `n_w` analytic functions of the wavelength and sweep `lam0` from 0.90 to
1.30 µm in steps of 0.02 µm:

    n_sio2(l) = sqrt(1 + 0.6961663*l^2/(l^2-0.0684043^2) + 0.4079426*l^2/(l^2-0.1162414^2) + 0.8974794*l^2/(l^2-9.896161^2))
    n_w(l)    = sqrt(1 + 0.75831*l^2/(l^2-0.01007) + 0.08495*l^2/(l^2-8.91377))

with `l` the wavelength in µm, and `n_a = n_w(l) + 0.016`. Export `n_eff(lam0)` and form
`beta2 = lam^3/(2*pi*c^2) * d2n_eff/dlam^2` by finite differences. Paper: **β₂ = −26.6 ps²/km** at 1.064 µm,
**β₃ = 0.135 ps³/km**, **D = +44.3 ps/(nm·km)**, zero-dispersion wavelength **0.927 µm**.

**5.3 Confinement loss against the number of rings.** Rebuild with the holes of rings 2–2 only (6 holes, one
ring) and rings 2–3 (18 holes, two rings), keeping the PML where it is. Paper: **10.5 dB/m** with one ring,
**2.6 × 10⁻⁵ dB/m** with two, below the numerical floor from three on. This is the run that shows five rings are
more than enough, and it is a good test of whether your PML is set up correctly: if the one-ring value is far
from 10 dB/m, the PML is reflecting.

**5.4 Index of the fundamental space-filling mode.** Separate small model: one rhombic unit cell of the lattice
(side Λ, angle 60°) with a single hole of diameter `dh` at its centre, periodic (Floquet) conditions with zero
wave vector on the four sides, no PML. The highest `n_eff` is `n_FSM`. Paper: **1.121** at 1.064 µm. This is the
number that proves the liquid core guides at all: `n_FSM` = 1.121 is well below the serum index 1.340, so the
core mode is bound.

## 6. If something goes wrong

| Symptom | Cause |
|---|---|
| `n_eff` comes out **above** `n_a` = 1.3399 | not a core mode. A mode guided in the liquid must have `n_eff < n_a`. Usually a PML or jacket mode; raise the number of modes and pick by core power fraction. |
| All modes near 1.44 | the solver found silica modes; check that the core material really overrides silica in the core domain. |
| Large imaginary part (10⁻³ or more) everywhere | a PML mode, or the PML is too thin or badly scaled. |
| `n_eff` differs from the paper in the third decimal | mesh in the 80 nm wall too coarse. Refine to 20–30 nm there and re-run. |
| Only one mode near 1.3320 | fine — the second polarization is degenerate; increase the number of modes to see it. |
| Core fraction around 90 % or less | the wall or a hole diameter is wrong; re-check `Rc + tw = 2.53 µm` and the six-tangent/six-truncated condition of Section 2.3. |

## 7. If COMSOL disagrees with the paper

Tell the other authors before the manuscript goes out. A disagreement is information, and it is much better to
find it now than in a referee report. The quantities to look at first, in order: the power fraction in the core
(geometry error), `n_FSM` (whether the cladding really has a low enough index), and only then `n_eff`.

If COMSOL agrees, you can add one sentence to Section 6.2 along the lines of

> The design was additionally checked with a commercial finite-element solver (COMSOL Multiphysics, Wave Optics
> Module), which gives n_eff = … and a core power fraction of … %, in agreement with the finite-difference
> values of Table 6.

and, if you want, the COMSOL mode plot as a second panel of Fig. 6. Report the numbers COMSOL actually gives.
