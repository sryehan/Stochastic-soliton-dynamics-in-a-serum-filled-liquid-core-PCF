% M2_EXACT_SOLUTIONS  Exact traveling-wave solutions (Section 4) and a numerical check that they satisfy
%   w^2 M'' - alpha1 M - 2 c M^3 = 0   (finite-difference residual).
clear; close all;
w = 1; k = 1.3; dth = 1e-3; th = (-6:dth:6)';
res = @(M, alpha1, c) max(abs(w^2*(M(3:end) - 2*M(2:end-1) + M(1:end-2))/dth^2 - alpha1*M(2:end-1) - 2*c*M(2:end-1).^3));
c = -0.7;  fprintf('sech : residual %.1e\n', res(k*w/sqrt(-c)*sech(k*th), k^2*w^2, c));
c = 0.7;   fprintf('tanh : residual %.1e\n', res(k*w/sqrt(c)*tanh(k*th), -2*k^2*w^2, c));
mu = 0.6;
[sn, cn, dn] = ellipj(k*th, mu);
c = -0.7;  fprintf('cn   : residual %.1e\n', res(sqrt(mu)*k*w/sqrt(-c)*cn, (2*mu - 1)*k^2*w^2, c));
           fprintf('dn   : residual %.1e\n', res(k*w/sqrt(-c)*dn, (2 - mu)*k^2*w^2, c));
c = 0.7;   fprintf('sn   : residual %.1e\n', res(sqrt(mu)*k*w/sqrt(c)*sn, -(1 + mu)*k^2*w^2, c));
% (the residuals are of the order of the discretization error dth^2 times the fourth derivative)

% rational (1/phi')-expansion: F = b/(r + a exp(-b theta)); family I reproduces the sech solution
b = 1.3; r = 0.8; a = 2.0; c = -0.7;
F = b./(r + a*exp(-b*th)); Q = b^2 - 2*b*r*F + 2*r^2*F.^2;
MI = b*w/sqrt(-c)*2*r*F.*(b - r*F)./Q; th0 = log(a/r)/b;
fprintf('family I minus sech form: %.1e\n', max(abs(MI - b*w/sqrt(-c)*sech(b*(th - th0)))));
c = 0.7; MIII = w/(2*sqrt(c))*(2*r*F - b);
fprintf('family III minus tanh form: %.1e\n', max(abs(MIII - b*w/(2*sqrt(c))*tanh(b*(th - th0)/2))));

% intensity maps of the bright and dark solitons (Fig. 4a,b)
x = linspace(-12, 12, 481); t = linspace(0, 10, 201); [X, T] = meshgrid(x, t);
Lb = bright_soliton(X, T, 0, 1, 0.5, 0.02, -1, 0);  Ld = dark_soliton(X, T, 0, 1, 0.5, 0.02, 1, 0);
figure('Position', [50 50 1000 380]);
subplot(1, 2, 1); imagesc(x, t, abs(Lb).^2); axis xy; colorbar; xlabel('x'); ylabel('t'); title('bright soliton |L|^2 (c = -1)');
subplot(1, 2, 2); imagesc(x, t, abs(Ld).^2); axis xy; colorbar; xlabel('x'); ylabel('t'); title('dark soliton |L|^2 (c = +1)');
print('-dpng', '-r150', 'm2_exact_solutions.png');
