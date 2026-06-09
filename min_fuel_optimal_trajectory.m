

clear; clc;

mu = 1;


B = [zeros(2,2); eye(2,2)];    

t0 = 0;
tf = 8;
tspan = [t0 tf];

u_max = 0.1;

rho_list = [1, 0.1, 1e-2, 1e-3];


lambda0_guess = [0.8; 0.1; 0.2; 1.1];


[Rm0, Vm0] = kep2cart(1.524, 0, pi, 0, 0, 0, mu);
x0_mars = [Rm0; Vm0];         


[Re0, Ve0] = kep2cart(1.0, 0, 0, 0, 0, 0, mu);
x0_earth = [Re0; Ve0];         % 4x1


opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',0.05);
[~, Xmars] = ode45(@(t,y) heliocentric_orbit(t,y,mu), tspan, x0_mars, opts);
x_mars_tf = Xmars(end,:).';    % 4x1

% fsolve options
fsolve_opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations',8000, ...
    'MaxIterations',200);

lambda0_sol_all = zeros(4, numel(rho_list));
exitflag_all = zeros(1, numel(rho_list));


for k_rho = 1:numel(rho_list)
    rho = rho_list(k_rho);

    fprintf('\n===== rho = %g =====\n', rho);

    Phi = @(lambda0) shooting_residual(lambda0, x0_earth, x_mars_tf, t0, tf, B, u_max, rho, mu);

    [lambda0_sol,~,exitflag,output] = fsolve(Phi, lambda0_guess, fsolve_opts);

    disp('---- fsolve status ----');
    disp(exitflag);
    disp(output);

    disp('---- lambda0 (solution) ----');
    disp(lambda0_sol);

    lambda0_sol_all(:,k_rho) = lambda0_sol;
    exitflag_all(k_rho) = exitflag;

    
    if exitflag > 0
        lambda0_guess = lambda0_sol;
    else
        
        fprintf('rho=%g did not converge; keeping previous lambda0_guess.\n', rho);
    end
end


rho = 1;
lambda0_sol = lambda0_sol_all(:,end);


z0 = [x0_earth; lambda0_sol];  % 8x1
[tt, ZZ] = ode45(@(t,z) state_dot(t,z,B,u_max,rho,mu), [t0 tf], z0, opts);

X   = ZZ(:,1:4);      
Lam = ZZ(:,5:8);       


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
        gamma = 0.5*u_max*(1 + tanh(S/rho));
        u = gamma*u_cap;
    end

    U(i,:) = u.';
    Gamma(i) = gamma;
    S_hist(i) = S;
end

x_tf = X(end,:).';
disp('---- terminal error x(tf)-x_mars(tf) ----');
disp(x_tf - x_mars_tf);


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

H_init = H_hist(1);
Hdiff_hist = H_hist - H_init;

figure;
plot(tt, Hdiff_hist, 'LineWidth', 1.5);
grid on;
legend('H(t) - H(0)');
xlabel('t'); ylabel('H(t) - H(0)');


figure;
plot(X(:,1), X(:,2), 'LineWidth', 1.5); hold on;
plot(Xmars(:,1), Xmars(:,2), '--', 'LineWidth', 1.5);
legend('Spacecraft (opt)','Mars (target)');
xlabel('x'); ylabel('y'); axis equal; grid on;

figure;
plot(tt, U(:,1), 'LineWidth', 1.5); hold on;
plot(tt, U(:,2), 'LineWidth', 1.5);
legend('u_x','u_y'); xlabel('t'); ylabel('u'); grid on;

figure;
plot(tt, Lam, 'LineWidth', 1.2);
xlabel('t'); ylabel('\lambda'); grid on;
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4');

figure;
plot(tt, Gamma, 'LineWidth', 1.5);
grid on;
xlabel('t'); ylabel('\Gamma');
title('\Gamma(t) (smoothed throttle)');

figure;
plot(tt, S_hist, 'LineWidth', 1.5);
grid on;
xlabel('t'); ylabel('S = ||B^T\lambda|| - 1');
title('Switching function');

t_ref = tt;        
x_ref = X;        
u_ref = U;        

save('ps4_reference_traj.mat', ...
     't_ref', 'x_ref', 'u_ref', ...
     'tf', 'mu', 'u_max');






function Phi = shooting_residual(lambda0, x0_earth, x_mars_tf, t0, tf, B, u_max, rho, mu)
    z0 = [x0_earth; lambda0];

    opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',0.05);

    try
        [~, Z] = ode45(@(t,z) state_dot(t,z,B,u_max,rho,mu), [t0 tf], z0, opts);

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

    a = -mu * r / r_mag^3;

    % Minimum-fuel (smoothed bang-bang) control
    p = B.'*lambda;            
    p_mag = norm(p);

    if p_mag < 1e-12
        u = zeros(2,1);
    else
        S = p_mag - 1;                              
        gamma = 0.5*u_max*(1 + tanh(S/rho));       
        u = gamma*u_cap;
    end

    x_dot = [v; a] + B*u;

   
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
