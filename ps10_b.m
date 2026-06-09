clear; clc; close all;

%% Load PS4 reference trajectory
load('ps4_reference_traj.mat');
% x_ref: Nref x 4 = [rx ry vx vy]
% u_ref: Nref x 2 = [ux uy]

%% Parameters
dt = 1e-2;
t_grid = 0:dt:tf;
N = length(t_grid);
M = 20;

sigma_pos  = 1e-2;
sigma_vel  = 1e-2;
sigma_dist = 1e-3;

P0 = diag([sigma_pos^2, sigma_pos^2, sigma_vel^2, sigma_vel^2]);

G = [0 0;
     0 0;
     sigma_dist 0;
     0 sigma_dist];

%% Interpolate reference trajectory/control onto PS10 grid
xbar = interp1(t_ref, x_ref, t_grid, 'linear');  % N x 4
ubar = interp1(t_ref, u_ref, t_grid, 'linear');  % N x 2

%% Monte Carlo simulation
Xmc = zeros(N, 4, M);

rng(1);
L0 = chol(P0, 'lower');

for m = 1:M

    Xmc(1,:,m) = xbar(1,:) + (L0*randn(4,1))';

    for k = 1:N-1

        xk = Xmc(k,:,m)';   % 4 x 1
        uk = ubar(k,:)';    % 2 x 1

        r = xk(1:2);
        v = xk(3:4);

        acc_grav = -mu*r/(norm(r)^3);

        f = [v;
             acc_grav + uk];

        dw = sqrt(dt)*randn(2,1);

        xnext = xk + f*dt + G*dw;

        Xmc(k+1,:,m) = xnext';
    end
end

%% Plot state histories
state_names = {'r_x','r_y','v_x','v_y'};

for i = 1:4
    figure; hold on; grid on;

    h_mc = plot(t_grid, squeeze(Xmc(:,i,1)), 'LineWidth', 0.8);

    for m = 2:M
        plot(t_grid, squeeze(Xmc(:,i,m)), 'LineWidth', 0.8);
    end

    h_ref = plot(t_grid, xbar(:,i), 'k', 'LineWidth', 2);

    xlabel('Time');
    ylabel(state_names{i});
    title(['Monte Carlo Time History of ', state_names{i}]);
    legend([h_mc h_ref], {'MC samples','Reference'});
end

%% Plot 2D position trajectories
figure; hold on; grid on; axis equal;

h_mc = plot(Xmc(:,1,1), Xmc(:,2,1), 'LineWidth', 0.8);

for m = 2:M
    plot(Xmc(:,1,m), Xmc(:,2,m), 'LineWidth', 0.8);
end

h_ref = plot(xbar(:,1), xbar(:,2), 'k', 'LineWidth', 2);

xlabel('r_x');
ylabel('r_y');
title('Monte Carlo Position Trajectories');
legend([h_mc h_ref], {'MC samples','Reference'});