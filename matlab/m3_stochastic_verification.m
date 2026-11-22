% M3_STOCHASTIC_VERIFICATION  Direct integration of the Stratonovich NLSE (no gauge transformation) and
% comparison with the exact stochastic soliton and with the exact statistics (Section 5):
%   |L|^2 unchanged,  |<L>| = exp(-varrho^2 t/2),  Var[arg L] = varrho^2 t.
clear; close all;
if exist('rng', 'file') || exist('rng', 'builtin'), rng(15197); else, randn('state', 15197); end
K = 1; m = 0; kappa = 0.02; c = -1;

% ---- (1) one realization, varrho = 0.3 ----------------------------------------------------------
varrho = 0.3; Lx = 60; N = 512; dt = 0.0025; Tend = 10; nt = round(Tend/dt);
x = (-N/2:N/2-1)*(Lx/N); k = 2*pi/Lx*[0:N/2-1, -N/2:-1];
u = bright_soliton(x, 0, 0, K, m, kappa, c, varrho); B = 0; err = 0; devI = 0; [dum, j0] = min(abs(x));
for it = 1:nt
    dB = sqrt(dt)*randn; B = B + dB;
    u = snlse_step_heun(u, k, dt, dB, kappa, c, varrho);
    ex = bright_soliton(x, it*dt, B, K, m, kappa, c, varrho);
    err = max(err, max(abs(u - ex))); devI = max(devI, abs(abs(u(j0))^2 - 1));
end
fprintf('single path, varrho = 0.3: max |L_num - L_exact| = %.1e, max | |L(0,t)|^2 - 1 | = %.1e\n', err, devI);

% ---- (2) ensemble statistics --------------------------------------------------------------------
np = 2000; Lx = 40; N = 256; dt = 0.01; nt = round(Tend/dt);
x = (-N/2:N/2-1)*(Lx/N); k = 2*pi/Lx*[0:N/2-1, -N/2:-1]; [dum, j0] = min(abs(x));
for varrho = [0.1 0.2 0.3]
    u = repmat(bright_soliton(x, 0, 0, K, m, kappa, c, varrho), np, 1); ph = zeros(np, 1); prev = zeros(np, 1);
    for it = 1:nt
        dB = sqrt(dt)*randn(np, 1);
        u = snlse_step_heun(u, k, dt, dB, kappa, c, varrho);
        a = angle(u(:, j0)*exp(-1i*0.5*(K^2 - m^2 - 2*kappa)*it*dt));       % phase relative to the noise-free solution
        ph = ph + angle(exp(1i*(a - prev))); prev = a;                        % unwrapped
    end
    n = 0.5*(K^2 - m^2 - 2*kappa);
    mf = abs(mean(u(:, j0)))/K*sqrt(-c);
    fprintf('varrho = %.1f: |<L>| = %.4f (theory %.4f), Var[arg L] = %.4f (theory %.4f +- %.4f), <|L|^2> = %.4f\n', ...
            varrho, mf, exp(-varrho^2*Tend/2), var(ph), varrho^2*Tend, varrho^2*Tend*sqrt(2/np), mean(abs(u(:, j0)).^2));
end

% ---- (3) the formal equation read in the Ito sense = Stratonovich equation with gain varrho^2/2 --
varrho = 0.3; u = bright_soliton(x, 0, 0, K, m, kappa, c, varrho); P0 = sum(abs(u).^2);
for it = 1:nt
    u = snlse_step_heun(u, k, dt, sqrt(dt)*randn, kappa, c, varrho, varrho^2/2);
end
fprintf('Ito reading, varrho = 0.3: power(t=10)/power(0) = %.4f, exp(varrho^2 t) = %.4f\n', sum(abs(u).^2)/P0, exp(varrho^2*Tend));
