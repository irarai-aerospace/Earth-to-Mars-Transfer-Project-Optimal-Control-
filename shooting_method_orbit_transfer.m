% =========================
% Problem Set 4 (a.2) – Single shooting, minimum-energy, fixed tf
% Planar (2D) heliocentric dynamics, nondimensional mu = 1
% Unknowns: lambda0 (4x1). Root condition: x(tf) - x_mars(tf) = 0
% =========================

clear; clc;

mu = 1;

% For planar state x = [rx; ry; vx; vy], control u is 2x1 acceleration
B = [zeros(2,2); eye(2,2)];    % 4x2

% initial guess for costate
for j = 1:20

lambda0_guess = unifrnd(-10,10,4,1);
% lambda0_guess = zeros(4,1);

t0 = 0;
tf = 8;
tspan = [t0 tf];

% -------- Planet initial conditions from given circular-orbit info --------
% NOTE: keplerian_to_cartesian must return PLANAR r(2x1), v(2x1) for these inputs.
% If your function returns 3D, modify it to output only x-y components.

% Mars: a=1.524, e=0, nu(t0)=pi
[Rm0, Vm0] = kep2cart(1.524, 0, pi, 0, 0, 0, mu);
x0_mars = [Rm0; Vm0];          % 4x1

% Earth: a=1.0, e=0, nu(t0)=0
[Re0, Ve0] = kep2cart(1.0, 0, 0, 0,  0, 0, mu);
x0_earth = [Re0; Ve0];         % 4x1

% Propagate Mars (uncontrolled) to get target x_mars(tf)
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~, Xmars] = ode45(@(t,y) heliocentric_orbit(t,y,mu), tspan, x0_mars, opts);
x_mars_tf = Xmars(end,:).';    % 4x1

% -------- Solve shooting problem for lambda0 --------
% Residual function: Phi(lambda0) = x(tf; lambda0) - x_mars(tf)
Phi = @(lambda0) shooting_residual(lambda0, x0_earth, x_mars_tf, t0, tf, B, mu);

% Solve for lambda0
fsolve_opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations',5000, ...
    'MaxIterations',200);

[lambda0_sol,~,exitflag,output] = fsolve(Phi, lambda0_guess, fsolve_opts);

if ~(exitflag > 0) || (output.iterations >= 50)
    fprintf('j=%d : not converged (exitflag=%d, iters=%d)\n', ...
            j, exitflag, output.iterations);
    continue;   % move to next random guess
end

fprintf('j=%d : converged (exitflag=%d, iters=%d)\n', ...
        j, exitflag, output.iterations);

% store solutions (use columns)
lambda0_sol_hist(:,j) = lambda0_sol;

disp('---- fsolve status ----');
disp(exitflag);
disp(output);

disp('---- lambda0 (solution) ----');
disp(lambda0_sol);

% -------- Integrate with solved lambda0 to report trajectories --------
z0 = [x0_earth; lambda0_sol];  % 8x1
[tt, ZZ] = ode45(@(t,z) state_dot(t,z,B,mu), [t0 tf], z0, opts);

X = ZZ(:,1:4);       % state history
Lam = ZZ(:,5:8);     % costate history

% control history u*(t) = -1/2 * B' * lambda(t)
U = zeros(length(tt),2);
for k = 1:length(tt)
    lambda_k = Lam(k,:).';
    U(k,:) = (-0.5*(B.'*lambda_k)).';
end
end
% Final state error check
x_tf = X(end,:).';
disp('---- terminal error x(tf)-x_mars(tf) ----');
disp(x_tf - x_mars_tf);
%%
% Hamiltonian evolution
v = x0_earth(3:4);
r = x0_earth(1:2);
r_mag = norm(r);

u = U(1,:).';                      % 2x1
f0 = [v; -mu*r/r_mag^3];           % 4x1
f  = f0 + B*u;                     % 4x1
H_init = (u.'*u) + lambda0_guess.'*f;    % scalar

H_hist     = zeros(length(tt),1);
Hdiff_hist = zeros(length(tt),1);

for i = 1:length(tt)
    r_t = X(i,1:2).'; 
    v_t = X(i,3:4).';
    r_t_mag = norm(r_t);

    u = U(i,:).';                                        % 2x1
    f_t = [v_t; -mu*r_t/r_t_mag^3] + B*u;                % 4x1

    lambda_t = Lam(i,:).';                               % 4x1

    H_hist(i)     = (u.'*u) + lambda_t.'*f_t;
    Hdiff_hist(i) = H_hist(i) - H_init;
end

figure;
plot(tt,Hdiff_hist,'LineWidth',1.5);
legend('H_{diff} = H(t) - H_{0}'); xlabel('t'); ylabel('H_{diff}')
%%
% (Optional) quick plots
figure; 

plot(X(:,1), X(:,2), 'LineWidth', 1.5); hold on;
plot(Xmars(:,1), Xmars(:,2), '--', 'LineWidth', 1.5);
legend('Spacecraft (opt)','Mars (target)'); xlabel('x'); ylabel('y'); axis equal; grid on;

figure; plot(tt, U(:,1), 'LineWidth', 1.5); hold on; plot(tt, U(:,2), 'LineWidth', 1.5);
legend('u_x','u_y'); xlabel('t'); ylabel('u'); grid on;

figure; plot(tt, Lam, 'LineWidth', 1.2);
xlabel('t'); ylabel('\lambda'); grid on; legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4');

% =========================
% Functions
% =========================

function Phi = shooting_residual(lambda0, x0_earth, x_mars_tf, t0, tf, B, mu)
    z0 = [x0_earth; lambda0];

    opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',0.05);

    try
        [~, Z] = ode45(@(t,z) state_dot(t,z,B,mu), [t0 tf], z0, opts);

        zf = Z(end,:).';
        xf = zf(1:4);

        Phi = xf - x_mars_tf;

        if any(~isfinite(Phi))
            Phi = 1e6*ones(4,1);
        end
    catch
        Phi = 1e6*ones(4,1);   % force fsolve to reject this guess
    end
end


function z_dot = state_dot(~, z, B, mu)
    % z = [x; lambda], where x=[r;v] in R^4, lambda in R^4
    x      = z(1:4);
    lambda = z(5:8);

    r = x(1:2);
    v = x(3:4);

    r_mag = norm(r);

    % Two-body accel (heliocentric): vdot = -mu r / r^3
    a = -mu * r / r_mag^3;

    % Optimal control for J = ∫ ||u||^2 dt  => u* = -1/2 B' lambda
    u = -0.5 * (B.' * lambda);

    x_dot = [v; a] + B*u;

    % Jacobian df0/dx for f0(x) = [v; a(r)]
    % da/dr = -mu( I/r^3 - 3 rr^T / r^5 )
    I2 = eye(2);
    da_dr = -mu * ( I2 / r_mag^3 - 3*(r*r.') / r_mag^5 );

    A = [zeros(2), I2;
         da_dr,    zeros(2)];

    % Costate dynamics: lambda_dot = -A' * lambda
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
%KEPLERIAN_TO_CARTESIAN  Planar/3D Keplerian elements to inertial Cartesian state.
% Inputs:
%   a   semi-major axis (nondimensional, e.g., AU units if l*=1 AU)
%   e   eccentricity
%   nu  true anomaly (rad)  <-- NOTE: for this assignment use nu_E=0, nu_M=pi
%   i   inclination (rad)
%   Om  RAAN (rad)
%   om  argument of periapsis (rad)
%   mu  gravitational parameter (nondimensional; mu=1 per assignment)
%
% Outputs:
%   R_inertial (2x1) position in inertial frame (planar x-y)
%   V_inertial (2x1) velocity in inertial frame (planar x-y)

    
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
