a0  = 9000e3;
mu  = 3.986e14;

Torb   = 2*pi*sqrt(a0^3/mu);

x0  = [a0; 0.2; deg2rad(50); deg2rad(10); deg2rad(20); 0];  % [a;e;i;sig;om;M]

% dt  = 1.0;                
% t0  = 0;
% tf  = 15*T;
% t   = (t0:dt:tf).';
% N   = numel(t);

% X   = zeros(N,6);
% X(1,:) = x0.';
% tk  = t0;
% yk  = x0;


tspan = [0, 15*Torb]; 

tf = 15*Torb;

RelTol = 1e-9;
AbsTol = 1e-9;
[T,X] = rk4_adaptive(@(t,x) gauss_variational_keplerian(t,x,mu), 0, tf, x0, 10, RelTol, AbsTol);

size(X)
length(T)

% for k = 2:N
%     xk = rk4singlestep(@(tt,xx) gauss_variational_keplerian(tt,xx,mu), dt, tk, yk);
%     tk = tk + dt;
%     yk = xk;              
%     X(k,:) = xk.';
% end

N = size(X,1);                 % MUST be this for adaptive
R1 = zeros(N,3);
V1 = zeros(N,3); 

a   = X(:,1);
e   = X(:,2);
inc = X(:,3);
sig = X(:,4);
om  = X(:,5);
M   = mod(X(:,6), 2*pi);  % wrap


for k = 1:N
    [Rk,Vk] = keplerian_to_cartesian(a(k), e(k), M(k), inc(k), sig(k), om(k), mu);
    R1(k,:) = Rk.'; 
    V1(k,:) = Vk.';
end

figure; plot3(R1(:,1),R1(:,2),R1(:,3)); axis equal; grid on

function dxdt = gauss_variational_keplerian(~,x,mu)
    a = x(1);
    n = sqrt(mu/a^3);
    dxdt = [0;0;0;0;0;n];
end
