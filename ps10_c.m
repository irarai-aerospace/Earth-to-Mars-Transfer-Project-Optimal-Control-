clear; clc; close all;

%% Load PS4 reference trajectory
load('ps4_reference_traj.mat');
% x_ref: Nref x 4
% u_ref: Nref x 2

%% Parameters
dt = 0.1;
t_grid = 0:dt:tf;
N = length(t_grid);

sigma_pos  = 1e-2;
sigma_vel  = 1e-2;
sigma_dist = 1e-3;

P0 = diag([sigma_pos^2, sigma_pos^2, sigma_vel^2, sigma_vel^2]);

G = [0 0;
     0 0;
     sigma_dist 0;
     0 sigma_dist];

opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Reference on covariance grid
xbar = interp1(t_ref, x_ref, t_grid, 'linear');  % N x 4
ubar = interp1(t_ref, u_ref, t_grid, 'linear');  % N x 2

%% Linear covariance propagation
P = zeros(4,4,N);
P(:,:,1) = P0;

Acell = cell(N-1,1);
Sigmacell = cell(N-1,1);

for k = 1:N-1

    xk = xbar(k,:)';
    uk = ubar(k,:)';

    % Continuous-time linearized A matrix
    Acont = continuousA(xk, mu);

    % Discrete transition matrix Phi
    Phi0 = eye(4);
    [~, PhiSol] = ode113(@(t,PhiVec) stmODE(t,PhiVec,Acont), ...
                         [t_grid(k), t_grid(k+1)], Phi0(:), opts_ode);

    Ak = reshape(PhiSol(end,:)',4,4);

    % Process-noise covariance Sigma_k from Eq. (5)
    Sigma0 = zeros(4,4);
    [~, SigmaSol] = ode113(@(t,SigmaVec) sigmaODE(t,SigmaVec,Acont,G), ...
                           [t_grid(k), t_grid(k+1)], Sigma0(:), opts_ode);

    Sigma_k = reshape(SigmaSol(end,:)',4,4);

    % Covariance update
    P(:,:,k+1) = Ak*P(:,:,k)*Ak' + Sigma_k;

    Acell{k} = Ak;
    Sigmacell{k} = Sigma_k;
end

%% Monte Carlo simulation on same dt = 0.1 grid
M = 20;
Xmc = zeros(N,4,M);
rng(1);

L0 = chol(P0,'lower');

for m = 1:M

    Xmc(1,:,m) = xbar(1,:) + (L0*randn(4,1))';

    for k = 1:N-1

        xk = Xmc(k,:,m)';
        uk = ubar(k,:)';

        % deterministic propagation with ZOH control
        [~, Xseg] = ode113(@(t,x) twoBodyEOM(x,uk,mu), ...
                           [t_grid(k), t_grid(k+1)], xk, opts_ode);

        xdet = Xseg(end,:)';

        % Brownian disturbance
        dw = sqrt(dt)*randn(2,1);

        xnext = xdet + G*dw;

        Xmc(k+1,:,m) = xnext';
    end
end

%% Plot deviations with 3-sigma bounds
state_names = {'r_x','r_y','v_x','v_y'};

for i = 1:4

    figure; hold on; grid on;

    % Monte Carlo deviations
    h_mc = plot(t_grid, squeeze(Xmc(:,i,1)) - xbar(:,i), ...
                'LineWidth', 0.8);

    for m = 2:M
        plot(t_grid, squeeze(Xmc(:,i,m)) - xbar(:,i), ...
             'LineWidth', 0.8);
    end

    % 3-sigma bound from covariance
    sigma_i = squeeze(sqrt(P(i,i,:)));

    h_sig1 = plot(t_grid,  3*sigma_i, 'k--', 'LineWidth', 2);
    h_sig2 = plot(t_grid, -3*sigma_i, 'k--', 'LineWidth', 2);

    xlabel('Time');
    ylabel(['\delta ', state_names{i}]);
    title(['Monte Carlo Deviation with 3\sigma Bound: ', state_names{i}]);

    legend([h_mc h_sig1], {'MC deviations','\pm3\sigma LinCov'});
end

%% ============================================================
%% Local functions

function xdot = twoBodyEOM(x,u,mu)

    r = x(1:2);
    v = x(3:4);

    xdot = [v;
            -mu*r/norm(r)^3 + u];
end

function A = continuousA(x,mu)

    r = x(1:2);
    rnorm = norm(r);

    dadr = -mu*(eye(2)/rnorm^3 - 3*(r*r')/rnorm^5);

    A = [zeros(2), eye(2);
         dadr,     zeros(2)];
end

function dPhiVec = stmODE(~,PhiVec,A)

    Phi = reshape(PhiVec,4,4);
    dPhi = A*Phi;

    dPhiVec = dPhi(:);
end

function dSigmaVec = sigmaODE(~,SigmaVec,A,G)

    Sigma = reshape(SigmaVec,4,4);

    dSigma = A*Sigma + Sigma*A' + G*G';

    dSigmaVec = dSigma(:);
end