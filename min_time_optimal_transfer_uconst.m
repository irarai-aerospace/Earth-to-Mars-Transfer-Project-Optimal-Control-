% =========================
% Problem Set 4 (b.2,b.3) – Single shooting, minimum-time, free tf
% Planar (2D) heliocentric dynamics, nondimensional mu = 1
% Unknowns: y = [lambda0; tf] (5x1)
% Root conditions: x(tf) - x_mars(tf) = 0 (4 eq) and H(tf)=0 (1 eq)
% Control constraint: ||u|| <= umax
% =========================

clear; clc;

mu = 1;

% Planar state x = [rx; ry; vx; vy], control u is 2x1 acceleration
B = [zeros(2,2); eye(2,2)];    % 4x2
umax = 0.1;

% initial guess (given)
lambda0_guess = [5;2;2;7];
t0 = 0;
tf_guess = 3;

% Unknown vector for fsolve
y_guess = [lambda0_guess; tf_guess];   % 5x1

% -------- Planet initial conditions (planar circular) --------
[Rm0, Vm0] = kep2cart(1.524, 0, pi, 0, 0, 0, mu);
x0_mars = [Rm0; Vm0];          % 4x1

[Re0, Ve0] = kep2cart(1.0, 0, 0, 0, 0, 0, mu);
x0_earth = [Re0; Ve0];         % 4x1

% -------- Solve shooting problem for y = [lambda0; tf] --------
Phi = @(y) shooting_residual(y, x0_earth, x0_mars, t0, B, umax, mu);

fsolve_opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations',8000, ...
    'MaxIterations',200);

[y_sol,~,exitflag,output] = fsolve(Phi, y_guess, fsolve_opts);

disp('---- fsolve status ----');
disp(exitflag);
disp(output);

lambda0_sol = y_sol(1:4);
tf_sol      = y_sol(5);

disp('---- lambda0 (solution) ----');
disp(lambda0_sol);

disp('---- tf (solution) ----');
disp(tf_sol);

% -------- Integrate with solved lambda0, tf to report trajectories --------
opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',0.05);

% Mars reference trajectory to tf_sol (for plotting/terminal check)
[tt_mars, Xmars] = ode45(@(t,y) heliocentric_orbit(t,y,mu), [t0 tf_sol], x0_mars, opts);
x_mars_tf = Xmars(end,:).';

% Spacecraft augmented trajectory
z0 = [x0_earth; lambda0_sol];  % 8x1
[tt, ZZ] = ode45(@(t,z) state_dot(t,z,B,mu,umax), [t0 tf_sol], z0, opts);

X   = ZZ(:,1:4);       % state history
Lam = ZZ(:,5:8);       % costate history

% control history (minimum-time, saturated)
U = zeros(length(tt),2);
for k = 1:length(tt)
    lambda_k = Lam(k,:).';
    p = B.'*lambda_k;                  % 2x1
    pmag = norm(p);
    if pmag < 1e-12
        u = zeros(2,1);
    else
        u = -umax*(p/pmag);
    end
    U(k,:) = u.';
end

% Final state error check
x_tf = X(end,:).';
disp('---- terminal error x(tf)-x_mars(tf) ----');
disp(x_tf - x_mars_tf);

% =========================
% (b.3) Hamiltonian evolution: plot H(t) - H(0)
% For minimum-time: H = 1 + lambda^T f(x,u)
% =========================
H_hist     = zeros(length(tt),1);
Hdiff_hist = zeros(length(tt),1);

for i = 1:length(tt)
    r_t = X(i,1:2).';
    v_t = X(i,3:4).';
    rmag = norm(r_t);

    u = U(i,:).';
    f_t = [v_t; -mu*r_t/rmag^3] + B*u;

    lambda_t = Lam(i,:).';
    H_hist(i) = 1 + lambda_t.'*f_t;
end

H_init = H_hist(1);
Hdiff_hist = H_hist - H_init;

figure;
plot(tt, Hdiff_hist, 'LineWidth', 1.5);
grid on;
legend('H(t) - H(0)');
xlabel('t'); ylabel('H(t) - H(0)');

% -------- Plots requested in (b.2) --------
figure;
plot(X(:,1), X(:,2), 'LineWidth', 1.5); hold on;
plot(Xmars(:,1), Xmars(:,2), '--', 'LineWidth', 1.5);
legend('Spacecraft (opt)','Mars (target)'); xlabel('x'); ylabel('y'); axis equal; grid on;

figure;
plot(tt, U(:,1), 'LineWidth', 1.5); hold on;
plot(tt, U(:,2), 'LineWidth', 1.5);
legend('u_x','u_y'); xlabel('t'); ylabel('u'); grid on;

figure;
plot(tt, Lam, 'LineWidth', 1.2);
xlabel('t'); ylabel('\lambda'); grid on;
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4');

% =========================
% Functions
% =========================

function Phi = shooting_residual(y, x0_earth, x0_mars, t0, B, umax, mu)
    % y = [lambda0; tf]
    lambda0 = y(1:4);
    tf      = y(5);

    % penalize non-physical tf
    if ~isfinite(tf) || tf <= 0
        Phi = 1e6*ones(5,1);
        return;
    end

    z0 = [x0_earth; lambda0];

    opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',0.05);

    try
        % Mars target at tf
        [~, Xmars] = ode45(@(t,ym) heliocentric_orbit(t,ym,mu), [t0 tf], x0_mars, opts);
        x_mars_tf = Xmars(end,:).';

        % Spacecraft augmented propagation to tf
        [~, Z] = ode45(@(t,z) state_dot(t,z,B,mu,umax), [t0 tf], z0, opts);

        zf = Z(end,:).';
        xf = zf(1:4);
        lambdaf = zf(5:8);

        % terminal matching
        Phi_x = xf - x_mars_tf;

        % H(tf) = 0 for minimum-time, free tf (autonomous, fixed terminal state)
        r = xf(1:2);
        v = xf(3:4);
        rmag = norm(r);

        % compute u(tf) from lambdaf
        p = B.'*lambdaf;
        pmag = norm(p);
        if pmag < 1e-12
            u = zeros(2,1);
        else
            u = -umax*(p/pmag);
        end

        f_tf = [v; -mu*r/rmag^3] + B*u;
        H_tf = 1 + lambdaf.'*f_tf;

        Phi = [Phi_x; H_tf];

        if any(~isfinite(Phi))
            Phi = 1e6*ones(5,1);
        end
    catch
        Phi = 1e6*ones(5,1);
    end
end

function z_dot = state_dot(~, z, B, mu, umax)
    % z = [x; lambda], x in R^4, lambda in R^4
    x      = z(1:4);
    lambda = z(5:8);

    r = x(1:2);
    v = x(3:4);

    r_mag = norm(r);

    % Uncontrolled dynamics
    a = -mu * r / r_mag^3;

    % Minimum-time optimal control with ||u||<=umax:
    % u* = -umax * (B'lambda)/||B'lambda||
    p = B.'*lambda;               % 2x1
    pmag = norm(p);
    if pmag < 1e-12
        u = zeros(2,1);
    else
        u = -umax * (p/pmag);
    end

    x_dot = [v; a] + B*u;

    % Jacobian df0/dx for f0(x) = [v; a(r)]
    I2 = eye(2);
    da_dr = -mu * ( I2 / r_mag^3 - 3*(r*r.') / r_mag^5 );

    A = [zeros(2), I2;
         da_dr,    zeros(2)];

    % Costate dynamics
    lambda_dot = -(A.') * lambda;

    z_dot = [x_dot; lambda_dot];
end

function dydt = heliocentric_orbit(~, y, mu)
    % Uncontrolled heliocentric two-body motion, planar 2D
    r = y(1:2);
    v = y(3:4);
    rmag = norm(r);

    drdt = v;
    dvdt = -mu * r / rmag^3;

    dydt = [drdt; dvdt];
end

function [R_inertial, V_inertial] = kep2cart(a,e,nu,i,Om,om,mu)
    if e == 0
        E = nu;
    else
        E = 2*atan2( sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2) );
    end

    rc = a*(1 - e*cos(E));

    r_pqw = rc * [cos(nu); sin(nu); 0];

    p = a*(1 - e^2);
    v_pqw = sqrt(mu/p) * [-sin(nu); (e + cos(nu)); 0];

    Rotation = rotz(rad2deg(Om)) * rotx(rad2deg(i)) * rotz(rad2deg(om));

    R3 = Rotation * r_pqw;
    V3 = Rotation * v_pqw;

    R_inertial = R3(1:2);
    V_inertial = V3(1:2);
end
