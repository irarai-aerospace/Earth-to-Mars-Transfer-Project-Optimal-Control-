a0  = 9000e3;
mu  = 3.986e14;

T   = 2*pi*sqrt(a0^3/mu);
[R0, V0] = keplerian_to_cartesian(a0,0.2,0,deg2rad(50),deg2rad(10),deg2rad(20),mu);
[p,f,g,h,k,L] = cartesian_to_equinoctial(R0,V0,mu);
x0  = [p;f;g;h;k;L];  

dt  = 5;                
t0  = 0;
tf  = 15*T;
t   = (t0:dt:tf).';
N   = numel(t);

X   = zeros(N,6);
X(1,:) = x0.';
tk  = t0;
yk  = x0;

tspan = [0, 15*T]; 
tf = 15*T;
RelTol = 1e-6;
AbsTol = 1e-6;
[T,X] = rk4_adaptive(@(t,x) gauss_variational_equinoctial(t,x,mu), 0, tf, x0, 10, RelTol, AbsTol);

% for j = 2:N
%     xk = rk4singlestep(@(tt,xx) gauss_variational_equinoctial(tt,xx,mu), dt, tk, yk);
%     tk = tk + dt;
%     yk = xk;              
%     X(j,:) = xk.';
% end
N = size(X,1);                 % MUST be this for adaptive
R = zeros(N,3);
V = zeros(N,3); 

p = X(:,1);
f = X(:,2);
g = X(:,3);
h = X(:,4);
k  = X(:,5);
L   = mod(X(:,6), 2*pi);  


for j = 1:N
    [Rk,Vk] = equinoctial_to_cartesian(p(j),f(j),g(j),h(j),k(j),L(j),mu);
    R(j,:) = Rk.'; 
    V(j,:) = Vk.';
end

figure; plot3(R(:,1),R(:,2),R(:,3)); axis equal; grid on

function dxdt = gauss_variational_equinoctial(~,x,mu)
    p = x(1);
    f = x(2);
    g = x(3);
    L = x(6);
    q = 1 + f*cos(L) + g*sin(L);
    dL = sqrt(mu/p^3) * q^2;
    dxdt = [0;0;0;0;0;dL];
end
