clc; clear; close all;

% Parameters
dt = 1e-3;
T  = 1;
N  = T/dt;              % number of increments
M  = 20;                % number of Monte Carlo samples
t  = (0:N-1)*dt;        % time associated with each increment

% Standard deviation of Brownian increment
sig_dw = sqrt(dt);
three_sigma = 3*sig_dw;

% Generate Brownian increments
% dw1 and dw2 are [M x N], each row = one sample path
dw1 = sig_dw * randn(M, N);
dw2 = sig_dw * randn(M, N);

% Plot component 1
figure;
plot(t, dw1', 'LineWidth', 0.8); hold on;
yline( three_sigma, 'r--', 'LineWidth', 1.5);
yline(-three_sigma, 'r--', 'LineWidth', 1.5);
xlabel('t');
ylabel('\Deltaw_{1,k}');
title('\Deltaw_{1,k} for M = 20 samples');
grid on;

% Plot component 2
figure;
plot(t, dw2', 'LineWidth', 0.8); hold on;
yline( three_sigma, 'r--', 'LineWidth', 1.5);
yline(-three_sigma, 'r--', 'LineWidth', 1.5);
xlabel('t');
ylabel('\Deltaw_{2,k}');
title('\Deltaw_{2,k} for M = 20 samples');
grid on;

% Optional: check empirical mean and std
fprintf('Component 1: empirical mean = %.5f, empirical std = %.5f\n', ...
    mean(dw1(:)), std(dw1(:)));
fprintf('Component 2: empirical mean = %.5f, empirical std = %.5f\n', ...
    mean(dw2(:)), std(dw2(:)));
fprintf('Theoretical std = %.5f, 3-sigma bound = %.5f\n', sig_dw, three_sigma);


%%
clc; clear; close all;

% Parameters
dt = 1e-3;
T  = 1;
N  = T/dt;                 % 1000 steps
M  = 20;                   % Monte Carlo samples
t  = 0:dt:T;               % state time vector, length N+1

sigma1 = 2;
sigma2 = 3;

G = [sigma1 0;
     0      sigma2];

% Preallocate state histories
x1 = zeros(M, N+1);
x2 = zeros(M, N+1);

% Monte Carlo simulation
for m = 1:M
    x = [0;0];             % initial state
    for k = 1:N
        dW = sqrt(dt) * randn(2,1);   % Brownian increment
        x = x + G*dW;                 % state update
        x1(m,k+1) = x(1);
        x2(m,k+1) = x(2);
    end
end

% Theoretical 3-sigma bounds
sig_x1 = sigma1 * sqrt(t);   % std of x1(t)
sig_x2 = sigma2 * sqrt(t);   % std of x2(t)

bound_x1 = 3 * sig_x1;
bound_x2 = 3 * sig_x2;

% Plot x1
figure;
plot(t, x1', 'LineWidth', 0.8); hold on;
plot(t,  bound_x1, 'r--', 'LineWidth', 2);
plot(t, -bound_x1, 'r--', 'LineWidth', 2);
xlabel('t');
ylabel('x_1(t)');
title('Monte Carlo histories of x_1(t)');
grid on;

% Plot x2
figure;
plot(t, x2', 'LineWidth', 0.8); hold on;
plot(t,  bound_x2, 'r--', 'LineWidth', 2);
plot(t, -bound_x2, 'r--', 'LineWidth', 2);
xlabel('t');
ylabel('x_2(t)');
title('Monte Carlo histories of x_2(t)');
grid on;

% Optional empirical check at final time
fprintf('At t = 1:\n');
fprintf('Empirical std of x1 = %.4f, Theoretical std = %.4f\n', std(x1(:,end)), sigma1);
fprintf('Empirical std of x2 = %.4f, Theoretical std = %.4f\n', std(x2(:,end)), sigma2);

