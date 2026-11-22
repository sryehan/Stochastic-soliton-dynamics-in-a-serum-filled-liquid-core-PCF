#!/bin/sh
# Regenerates all results and figures.  Run from the repository root:  sh run_all.sh
set -e
cd code
python verify_math.py      > ../out/verify_math.log
python verify_families.py  > ../out/verify_families.log
python sim_pathwise.py
python sim_final.py path
python sim_final.py ens
python sim_final.py ito
python validate_he11.py     > ../out/validate_he11.log
python validate_six_hole.py > ../out/validate_six_hole.log
for stage in design sweep sens conv modes rings fsm; do
    python fibre_results.py $stage > ../out/run_$stage.log
done
# glucose data set: download the CSV first (see data/README.md)
if [ -f ../data/diabetes_prediction_dataset.csv ]; then
    python kaggle_summary.py ../data/diabetes_prediction_dataset.csv
fi
python fibre_results.py kaggle > ../out/run_kaggle.log
python derive.py
python sim_fibre.py > ../out/run_prop.log
python kaggle_analysis.py
python validate_bpm_scalar.py > ../out/validate_bpm_scalar.log
# beam propagation (Section 6.5): run matlab/m6_beam_propagation.m in MATLAB or Octave first
if [ -f ../matlab/bpm_results.mat ]; then python bpm_export.py ../matlab/bpm_results.mat; fi
python bpm_movie.py            # videos/ and out/bpm_side_view.npz (needs ffmpeg)
for f in fig_phase fig_solutions fig_stoch fig_fibre fig_bpm fig_side fig_sensing fig_kaggle; do python $f.py; done
python export_figure_data.py   # out/figure_data.mat for matlab/paper_figures.m
python gen_matlab.py
