# MATLAB scripts

MATLAB is the main environment of this repository for drawing and editing the figures. The scripts use only base
MATLAB (R2016b or later) and also run in GNU Octave.

## Quick start: the figures of the paper as editable MATLAB figures

    cd matlab
    paper_figures

`paper_figures.m` reads `../out/figure_data.mat` (the computed results of the paper, exported by
`code/export_figure_data.py`) and draws Figs. 5–12. Each figure is saved in `matlab/figures_paper/` as a `.fig` file,
which you can open in MATLAB and edit (colours, labels, fonts, axis limits, line widths), and as a 300 dpi `.png`.
It takes a few seconds and needs no other toolbox. Figs. 2–4 (phase portraits, bifurcation, exact solutions) are drawn
by `m1_phase_portraits.m` and `m2_exact_solutions.m`; Fig. 1 (geometry and set-up sketch) exists only in Python
(`code/fig_fibre.py`).

## All scripts

| Script | Content | Paper |
|---|---|---|
| `paper_figures.m` | draws Figs. 5–12 from the computed data; saves editable `.fig` and `.png` in `figures_paper/` (seconds) | Figs. 5–12 |
| `m1_phase_portraits.m` | classification of the equilibria and the six phase portraits | Section 3, Fig. 2 |
| `m2_exact_solutions.m` | exact solutions, finite-difference check of the ODE, families of the (1/phi')-expansion | Section 4, Fig. 4 |
| `m3_stochastic_verification.m` | direct Stratonovich simulation (split-step + stochastic Heun), one path, 2000-path ensemble, Ito reading | Section 5, Fig. 5, Table 4 |
| `m4_fiber_propagation.m` | soliton with absorption and third-order dispersion against linear propagation | Section 6.4, Fig. 8 |
| `m5_sensing_iso15197.m` | sensitivities, detection limit, ISO 15197:2013 test on the Kaggle glucose levels | Section 7, Figs. 11, 12 |
| `m6_beam_propagation.m` | beam propagation along three liquid-core fibers (liquid core in silica, conventional cladding, this work), scalar mode of this work, phase difference between two glucose levels | Section 6.5, Fig. 9 |
| `m7_propagation_movie.m` | the same three fibers advanced together; writes 4 videos (A, B, C, and A \| B \| C) and the side view to `videos_m7/` (MPEG-4, 30 fps, up to 1920 x 1080; 6–10 min) | Fig. 10, Supplementary Videos 1–4 |

Helper functions: `snlse_step_heun.m`, `bright_soliton.m`, `dark_soliton.m`, `pulse_fwhm.m`, `ifelse_str.m`;
for `m6` and `m7`: `pcf_index_map.m` (cross-sections), `bpm_run.m` (split-step propagation), `bpm_mode.m` (mode by imaginary-distance
propagation), `capillary_leaky_mode.m` (exact leaky mode of a liquid core in silica).
`m4` and `m5` are written by `code/gen_matlab.py`, which inserts the modal parameters from `out/numbers.json`.

The deterministic scripts m1–m5 reproduce the Python numbers: for example `m4` gives pulse-width ratios 1.027, 1.204 and
1.484 at 2, 5 and 8 cm for the 30 fs soliton, and `m5` gives a temperature bound of 0.135 K. The stochastic script
uses different random numbers from the Python run, so it agrees statistically, not digit by digit. The Octave run of
m1–m5 recorded in `octave_output.txt` was headless and shows the text output only. The full-vector mode solver of
Section 6 exists only in Python (`code/pcf7.py`, `code/fibre_results.py`); its results are in `out/` and are what
`paper_figures.m` plots.

## m6: beam propagation (Section 6.5)

`m6_beam_propagation.m` solves the paraxial wave equation (Eq. 34 of the paper)

    d(psi)/dz = i/(2 k0 n0) (d2/dx2 + d2/dy2) psi + i k0 (n(x,y)^2 - n0^2)/(2 n0) psi

over the full cross-section with the split-step Fourier method; the distance z plays the role of time. It was run in
GNU Octave 8.4 with the gnuplot toolkit (9 min on two cores); the text output is in `octave_output_m6.txt` and the four
figures it writes are in `figures_m6/`. It also saves `bpm_results.mat` (25 MB, not stored in the repository).
`code/bpm_export.py` and `code/fig_bpm.py` turn that file into Fig. 9 of the paper (Python), so a MATLAB run can replace the figure:

    cd matlab;  run m6_beam_propagation          % in MATLAB or Octave
    cd ../code; python bpm_export.py ../matlab/bpm_results.mat; python fig_bpm.py

Values obtained (grid 640 x 640, pixel 50 nm, step 25 nm):

| Quantity | Value |
|---|---|
| A, liquid core in silica: decay of the core power | 0.0293 per µm = 127 dB/mm (exact leaky mode 0.0294 per µm) |
| B, conventional cladding (d/Λ = 0.5): core power at 50 / 200 µm | 0.144 / 0.023 (cladding index n_FSM = 1.440 > n_a) |
| C, this work: core power from 100 to 200 µm | 0.974 ± 0.004; 0.990 of the power still in the window at 200 µm |
| C, scalar mode | n_eff = 1.332434 (independent finite-difference check 1.332442; full-vector 1.332018), 98.8 % in the liquid, A_eff = 12.41 µm² |
| C, 300 against 80 mg/dL | 1.808 mrad/µm, S_n = 0.994 (full-vector 0.9964), 90.4 rad over 5 cm (paper 90.6 rad) |

The model is scalar and paraxial. Set `fast = true` to halve the run time (step 50 nm, about 1 % numerical loss per
100 µm). The step must stay small because of the silica–air index jumps.

## COMSOL

The COMSOL model moved to the `comsol/` folder of this repository: `lcpcf_model.m` (LiveLink for MATLAB),
`COMSOL_guide.md` (step-by-step build in the COMSOL Desktop, with the values to compare against) and the geometry
files. None of it has been executed by the authors of this repository; see that folder's guide.
