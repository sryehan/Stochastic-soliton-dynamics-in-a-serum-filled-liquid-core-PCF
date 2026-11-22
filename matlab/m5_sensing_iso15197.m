% M5_SENSING_ISO15197  Phase and group-delay sensitivity to glucose, noise-limited detection limit and the virtual
% accuracy test against ISO 15197:2013 on the glucose levels of the Kaggle "Diabetes prediction dataset"
% (Section 7).  The level counts below are the aggregate of the column blood_glucose_level (out/kaggle_levels.json).
clear; close all;
lam0 = 1.064e-6; k0 = 2*pi/lam0; c0 = 299792458;
Sn = 0.996369; dngdna = 0.996733;                 % d n_eff/d n_a and d n_g/d n_a from the mode solver
LD = 0.033824; Lf = 0.05;                      % dispersion length for T0 = 30 fs (m), fiber length (m)
dndC = 1.4e-06; dndT = 0.0001;              % RIU per mg/dL (glucose), RIU per K (water, magnitude)
levels = [80 85 90 100 126 130 140 145 155 158 159 160 200 220 240 260 280 300];
counts = [7106 6901 7112 7025 7702 7794 7732 7679 7575 7026 7759 7712 7600 603 636 635 729 674];
SphiC = k0*Lf*Sn*dndC; StauC = Lf/c0*dngdna*dndC*1e15;
fprintf('phase sensitivity %.3e rad/RIU = %.4f rad per mg/dL; delay sensitivity %.1f ps/RIU = %.4f fs per mg/dL\n', k0*Lf*Sn, SphiC, Lf/c0*dngdna*1e12, StauC);
fprintf('one fringe = %.2f mg/dL = %.2f fs of group delay; temperature cross-sensitivity %.1f mg/dL per K\n', 2*pi/SphiC, 2*pi/SphiC*StauC, dndT/dndC);
for varrho = [0.01 0.1 0.3]
    fprintf('varrho = %.2f: 3-sigma detection limit %.3f mg/dL\n', varrho, 3*varrho*sqrt(Lf/LD)/SphiC);
end
w = counts/sum(counts); tol = 15*(levels < 100) + 0.15*levels.*(levels >= 100);
fprintf('records %d, mean %.2f, SD %.2f mg/dL, share below 100 mg/dL %.3f\n', sum(counts), sum(w.*levels), sqrt(sum(w.*(levels - sum(w.*levels)).^2)), sum(w(levels < 100)));
frac = @(s) sum(w.*erf(tol/(sqrt(2)*s)));                                   % fraction of records within the ISO limits
sC = @(varrho, sT) sqrt((varrho*sqrt(Lf/LD)/SphiC)^2 + (dndT/dndC*sT)^2);   % standard deviation of the reading (mg/dL)
for sT = [0.05 0.10 0.15 0.20]
    fprintf('varrho = 0.1, sigma_T = %.2f K: %.2f %% of the records within the ISO 15197:2013 limits\n', sT, 100*frac(sC(0.1, sT)));
end
s95 = fzero(@(s) frac(s) - 0.95, [1 50]);
fprintf('95 %% requirement: total SD <= %.2f mg/dL, i.e. sigma_T <= %.4f K (varrho = 0.3), %.4f K (varrho = 1)\n', s95, ...
        sqrt(s95^2 - (0.3*sqrt(Lf/LD)/SphiC)^2)/(dndT/dndC), sqrt(s95^2 - (1.0*sqrt(Lf/LD)/SphiC)^2)/(dndT/dndC));
