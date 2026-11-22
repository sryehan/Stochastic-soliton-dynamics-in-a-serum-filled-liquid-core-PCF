% M1_PHASE_PORTRAITS  Phase portraits of  M' = P,  P' = (alpha1*M + beta1*M^3)/w^2  (Section 3, Fig. 2).
% alpha1 = m^2 + 2(n + kappa), beta1 = 2c.  The orbits are level curves of the Hamiltonian H(M,P).
clear; close all;
w = 1;
cases = [0.5 -0.8 0.02 -0.5;  1 -0.5 0 -0.5;  1 0.3 0.1 -1;  0.3 -1.5 0.01 2;  1 -0.5 0 0.5;  1.5 0.5 0.12 0.5];   % m, n, kappa, c
names = {'(a) global center', '(b) nonlinear center', '(c) saddle + two centers: bright soliton', ...
         '(d) center + two saddles: dark soliton', '(e) degenerate saddle', '(f) saddle only'};
% ---- classification of the equilibria (eigenvalues +-sqrt(alpha1 + 3*beta1*Me^2)/w) ------------------------
for j = 1:6
    m = cases(j,1); n = cases(j,2); kappa = cases(j,3); c = cases(j,4);
    alpha1 = m^2 + 2*(n + kappa); beta1 = 2*c;
    fprintf('case %d: alpha1 = %6.2f, beta1 = %5.1f, origin: %s', j, alpha1, beta1, ifelse_str(alpha1 > 0, 'saddle', alpha1 < 0, 'center', 'degenerate'));
    if alpha1*beta1 < 0
        fprintf('; E+- = (+-%.4f, 0): %s', sqrt(-alpha1/beta1), ifelse_str(alpha1 > 0, 'centers', true, 'saddles', ''));
    end
    fprintf('\n');
end
% ---- phase portraits ---------------------------------------------------------------------------------------
[M, P] = meshgrid(linspace(-2.2, 2.2, 401));
figure('Position', [50 50 1200 760]);
for j = 1:6
    m = cases(j,1); n = cases(j,2); kappa = cases(j,3); c = cases(j,4);
    alpha1 = m^2 + 2*(n + kappa); beta1 = 2*c;
    H = P.^2/2 - (alpha1*M.^2/2 + beta1*M.^4/4)/w^2;
    subplot(2, 3, j); contour(M, P, H, 40); hold on;
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
    if alpha1*beta1 < 0
        Me = sqrt(-alpha1/beta1); plot([-Me Me], [0 0], 'ks', 'MarkerFaceColor', 'r');
        He = -(alpha1*Me^2/2 + beta1*Me^4/4)/w^2;
        if alpha1 > 0, contour(M, P, H, [0 0], 'b', 'LineWidth', 2); else, contour(M, P, H, [He He], 'r', 'LineWidth', 2); end
    end
    xlabel('M'); ylabel('P = dM/d\theta'); title(sprintf('%s\n\\alpha_1 = %.2f, \\beta_1 = %.1f', names{j}, alpha1, beta1));
end
print('-dpng', '-r150', 'm1_phase_portraits.png');
