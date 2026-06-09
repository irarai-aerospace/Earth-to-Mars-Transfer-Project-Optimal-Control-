% =========================
% Optional (extra credit) – Random initial guesses for minimum-fuel (fixed tf)
% Perform N=20 random starts, each element of lambda0 ~ U(-10,10)
% Use the SAME continuation-on-rho strategy for each random start.
% Discuss convergence sensitivity (expected higher than min-energy/min-time).
%
% Output:
%   - converged vs not
%   - iterations / exitflags / residual norms
%   - (optional) clusters of lambda0 solutions
% =========================

clear; clc; close all;

mu = 1;

% Planar state x = [rx; ry; vx; vy], control u is 2x1 acceleration
B = [zeros(2,2); eye(2,2)];    % 4x2

t0 = 0;
tf = 8;
tspan = [t0 tf];

u_max = 0.1;

% Continuation list (as given)
rho_list = [1, 0.1, 1e-2, 1e-3];
rho_final = rho_list(end);

% -------- Planet initial conditions (planar circular) --------
[Rm0, Vm0] = kep2cart(1.524, 0, pi, 0, 0, 0, mu);
x0_mars = [Rm0; Vm0];

[Re0, Ve0] = kep2cart(1.0, 0, 0, 0, 0, 0, mu);
x0_earth = [Re0; Ve0];

% Propagate Mars (uncontrolled) to get target x_mars(tf)
opts_orbit = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.05);
[~, Xmars] = ode45(@(t,y) heliocentric_orbit(t,y,mu), tspan, x0_mars, opts_orbit);
x_mars_tf = Xmars(end,:).';

% fsolve options (assignment tolerance 1e-8)
fsolve_opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-8, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations',8000, ...
    'MaxIterations',200);

% =========================
% Random-start study
% =========================
N = 20;

lambda0_guess_hist = nan(4,N);
lambda0_sol_hist   = nan(4,N);
exitflag_hist      = nan(1,N);
iters_total_hist   = nan(1,N);      % total fsolve iterations across rho continuation
resnorm_hist       = nan(1,N);      % terminal residual norm at final rho
rho_fail_hist      = nan(1,N);      % rho at which it failed (if failed)

fprintf('===== Minimum-fuel random-start study (N=%d) =====\n', N);

for j = 1:N
    % Random initial guess per extra-credit statement
    lambda0_guess = unifrnd(-10,10,4,1);
    lambda0_guess_hist(:,j) = lambda0_guess;

    lambda_guess = lambda0_guess;
    total_iters = 0;
    exitflag_last = -Inf;
    failed_at_rho = NaN;

    % Continuation over rho for THIS random start
    for k_rho = 1:numel(rho_list)
        rho = rho_list(k_rho);

        Phi = @(lambda0) shooting_residual(lambda0, x0_earth, x_mars_tf, t0, tf, B, u_max, rho, mu);

        try
            [lambda_sol,~,exitflag,output] = fsolve(Phi, lambda_guess, fsolve_opts);
            total_iters = total_iters + output.iterations;

            % Check residual too (sometimes exitflag>0 but residual still not small)
            res = Phi(lambda_sol);
            resnorm = norm(res);

            if ~(exitflag > 0) || ~isfinite(resnorm) || resnorm > 1e-3
                exitflag_last = exitflag;
                failed_at_rho = rho;
                break;
            end

            % Accept solution and use as next rho guess
            lambda_guess = lambda_sol;
            exitflag_last = exitflag;

        catch
            exitflag_last = -999;
            failed_at_rho = rho;
            break;
        end
    end

    % Store results
    iters_total_hist(j) = total_iters;
    exitflag_hist(j)    = exitflag_last;
    rho_fail_hist(j)    = failed_at_rho;

    if isfinite(failed_at_rho)
        % failed during continuation
        lambda0_sol_hist(:,j) = nan(4,1);
        resnorm_hist(j) = inf;
        fprintf('run %2d/%2d: FAIL at rho=%g (exitflag=%d, total iters=%d)\n', ...
            j, N, failed_at_rho, exitflag_last, total_iters);
    else
        % success through all rhos
        lambda0_sol_hist(:,j) = lambda_guess;

        % final residual at rho_final
        res = shooting_residual(lambda_guess, x0_earth, x_mars_tf, t0, tf, B, u_max, rho_final, mu);
        resnorm_hist(j) = norm(res);

        fprintf('run %2d/%2d: OK   (total iters=%d, ||Phi||=%.3e)\n', ...
            j, N, total_iters, resnorm_hist(j));
    end
end

% Convergence definition (final residual small)
tol_res = 1e-6;
converged = isfinite(resnorm_hist) & (resnorm_hist < tol_res);

fprintf('\n========== SUMMARY (Minimum-fuel random starts) ==========\n');
fprintf('Converged: %d / %d (tol ||Phi|| < %.1e)\n', sum(converged), N, tol_res);

if any(converged)
    fprintf('Total iterations (converged): min=%d, median=%d, max=%d\n', ...
        min(iters_total_hist(converged)), round(median(iters_total_hist(converged))), max(iters_total_hist(converged)));
    fprintf('Residual norms (converged): min=%.2e, median=%.2e, max=%.2e\n', ...
        min(resnorm_hist(converged)), median(resnorm_hist(converged)), max(resnorm_hist(converged)));
end

% =========================
% Plots for discussion
% =========================
figure;
stem(1:N, converged, 'filled'); ylim([-0.1 1.1]);
xlabel('random start index'); ylabel('converged (1=yes,0=no)');
grid on; title('Minimum-fuel: convergence outcomes (random \lambda_0)');

figure;
semilogy(1:N, resnorm_hist, 'o', 'LineWidth', 1.2); grid on; hold on;
yline(tol_res,'--','tol');
xlabel('random start index'); ylabel('||\Phi(\lambda_0)|| at final \rho');
title('Minimum-fuel: terminal residual norms across random starts');

figure;
plot(1:N, iters_total_hist, 'o', 'LineWidth', 1.2); grid on;
xlabel('random start index'); ylabel('total fsolve iterations (sum over \rho)');
title('Minimum-fuel: total iterations across random starts');

% Optional: show where failures happen in continuation
figure;
rho_fail_plot = rho_fail_hist;
rho_fail_plot(isnan(rho_fail_plot)) = 0; % 0 means success all the way
stem(1:N, rho_fail_plot, 'filled');
grid on;
xlabel('random start index'); ylabel('\rho where it failed (0 = success)');
title('Minimum-fuel: continuation failure location');

% =========================
% If you want ONE representative successful trajectory plot:
% Pick the first converged run and integrate full history at rho_final.
% =========================
idx_ok = find(converged, 1, 'first');
if ~isempty(idx_ok)
    lambda0_sol = lambda0_sol_hist(:,idx_ok);

    opts_state = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.05);
    z0 = [x0_earth; lambda0_sol];
    [tt, ZZ] = ode45(@(t,z) state_dot(t,z,B,u_max,rho_final,mu), [t0 tf], z0, opts_state);

    X = ZZ(:,1:4);
    Lam = ZZ(:,5:8);

    % Control history
    U = zeros(length(tt),2);
    Gamma = zeros(length(tt),1);
    S_hist = zeros(length(tt),1);

    for i = 1:length(tt)
        lambda_t = Lam(i,:).';
        p = B.'*lambda_t;
        p_mag = norm(p);

        if p_mag < 1e-12
            u = zeros(2,1);
            gamma = 0;
            S = -1;
        else
            u_cap = -p/p_mag;
            S = p_mag - 1;
            gamma = 0.5*u_max*(1 + tanh(S/rho_final));
            u = gamma*u_cap;
        end

        U(i,:) = u.';
        Gamma(i) = gamma;
        S_hist(i) = S;
    end

    % Hamiltonian evolution H(t)-H(0): H = Gamma + lambda^T (f0 + Bu)
    H_hist = zeros(length(tt),1);
    for i = 1:length(tt)
        r_t = X(i,1:2).';
        v_t = X(i,3:4).';
        rmag = norm(r_t);

        u = U(i,:).';
        f_t = [v_t; -mu*r_t/rmag^3] + B*u;

        lambda_t = Lam(i,:).';
        H_hist(i) = Gamma(i) + lambda_t.'*f_t;
    end

    figure;
    plot(tt, H_hist - H_hist(1), 'LineWidth', 1.5);
    grid on; xlabel('t'); ylabel('H(t)-H(0)');
    title(sprintf('Minimum-fuel: H(t)-H(0) (example run %d)', idx_ok));

    figure;
    plot(X(:,1), X(:,2), 'LineWidth', 1.5); hold on;
    plot(Xmars(:,1), Xmars(:,2), '--', 'LineWidth', 1.5);
    legend('Spacecraft (opt)','Mars (target)');
    xlabel('x'); ylabel('y'); axis equal; grid on;
    title(sprintf('Minimum-fuel: trajectory (example run %d)', idx_ok));

    figure;
    plot(tt, Gamma, 'LineWidth', 1.5); grid on;
    xlabel('t'); ylabel('\Gamma'); title('Smoothed throttle \Gamma(t)');

    figure;
    plot(tt, S_hist, 'LineWidth', 1.5); grid on;
    xlabel('t'); ylabel('S = ||B^T\lambda|| - 1'); title('Switching function S(t)');
end

% =========================
% Functions (same dynamics as your c.3/c.4 code)
% =========================

function Phi = shooting_residual(lambda0, x0_earth, x_mars_tf, t0, tf, B, u_max, rho, mu)
    z0 = [x0_earth; lambda0];

    % Event: stop obviously bad trajectories to avoid hanging on random guesses
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',0.05, ...
                  'Events', @(t,z) stop_bad_traj(t,z));

    try
        [~, Z, te] = ode45(@(t,z) state_dot(t,z,B,u_max,rho,mu), [t0 tf], z0, opts);

        if ~isempty(te) || isempty(Z)
            Phi = 1e6*ones(4,1);
            return;
        end

        zf = Z(end,:).';
        xf = zf(1:4);

        Phi = xf - x_mars_tf;

        if any(~isfinite(Phi))
            Phi = 1e6*ones(4,1);
        end
    catch
        Phi = 1e6*ones(4,1);
    end
end

function z_dot = state_dot(~, z, B, u_max, rho, mu)
    x      = z(1:4);
    lambda = z(5:8);

    r = x(1:2);
    v = x(3:4);
    r_mag = norm(r);

    % Uncontrolled accel
    a = -mu * r / r_mag^3;

    % Minimum-fuel (smoothed bang-bang) control
    p = B.'*lambda;
    p_mag = norm(p);

    if p_mag < 1e-12
        u = zeros(2,1);
    else
        u_cap = -p/p_mag;
        S = p_mag - 1;
        gamma = 0.5*u_max*(1 + tanh(S/rho));
        u = gamma*u_cap;
    end

    x_dot = [v; a] + B*u;

    % Costate dynamics
    I2 = eye(2);
    da_dr = -mu * ( I2 / r_mag^3 - 3*(r*r.') / r_mag^5 );

    A = [zeros(2), I2;
         da_dr,    zeros(2)];

    lambda_dot = -(A.') * lambda;

    z_dot = [x_dot; lambda_dot];
end

function dydt = heliocentric_orbit(~, y, mu)
    r = y(1:2);
    v = y(3:4);
    rmag = norm(r);

    dydt = [v; -mu * r / rmag^3];
end

function [value, isterminal, direction] = stop_bad_traj(~, z)
    x = z(1:4);
    r = x(1:2);
    rmag = norm(r);

    rmin = 0.2;
    rmax = 10.0;

    bad = (~all(isfinite(z))) || (rmag < rmin) || (rmag > rmax);

    value = double(~bad);
    isterminal = 1;
    direction = 0;
end

function [R_inertial, V_inertial] = kep2cart(a,e,nu,i,Om,om,mu)
% Planar Keplerian -> Cartesian (x-y only); signature unchanged.
    %#ok<INUSD>
    p = a*(1 - e^2);
    r = p/(1 + e*cos(nu));

    R_inertial = r*[cos(nu); sin(nu)];
    V_inertial = sqrt(mu/p)*[-sin(nu); e + cos(nu)];
end