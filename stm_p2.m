% ====================== 2BP + J2 STM integration (15 orbits) ======================
clear; clc;

% ---- constants (SI) ----
mu = 3.986e14;          % m^3/s^2
J2 = 1.08263e-3;        % (use standard value)
Re = 6378.1363e3;       % m

a  = 9000e3;            % m
T  = 2*pi*sqrt(a^3/mu); % s
t0 = 0;
tf = 15*T;
a_pert = a + 9000e3 + 1e3;
omega_pert = deg2rad(10) + deg2rad(10);
% ---- initial condition from Kepler elements (example you had) ----
% keplerian_to_cartesian(a,e,Omega,i,omega,M,mu)  (adjust if your signature differs)
[R,V] = keplerian_to_cartesian(9000e3, 0.2, 0, deg2rad(50), deg2rad(10), deg2rad(20), mu);

x0 = [R(:); V(:)];      % 6x1

% ---- build augmented initial condition with STM = I ----
n  = 6;
Phi0 = eye(n);
X0_aug = [x0; Phi0(:)];

% ---- build numeric dynamics f(x) and A(x)=df/dx from symbolic ONCE ----
syms x y z vx vy vz real
r  = [x; y; z];
v  = [vx; vy; vz];
r2 = x^2 + y^2 + z^2;
rm = sqrt(r2);

% 2-body accel
a2b = -mu * r / rm^3;

% J2 accel (ECI, using z along Earth's spin axis)
fac = (3/2)*J2*mu*Re^2 / rm^5;
z2r2 = (z^2)/r2;
aJ2 = fac * [ x*(5*z2r2 - 1);
             y*(5*z2r2 - 1);
             z*(5*z2r2 - 3) ];

a = a2b + aJ2;

f_sym = [ v;
          a ];

vars = [x y z vx vy vz];
A_sym = jacobian(f_sym, vars);

% numeric function handles (NO time dependence here)
f_fun = matlabFunction(f_sym, 'Vars', {x,y,z,vx,vy,vz});
A_fun = matlabFunction(A_sym, 'Vars', {x,y,z,vx,vy,vz});


% ---- integrate augmented system ----
N = 2000;   
this = linspace(t0,tf,N); % only need final STM; change to linspace(t0,tf,N) if you want history
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

[traj,Phi] = STM_compute(t0, tf, this, X0_aug, n, f_fun, A_fun, opts);

Xhis = deval(traj, this);        % (6+36) x N
xhis = Xhis(1:n, :);             % 6 x N

% ---- xdot history at the same epochs ----
xdot = zeros(n, numel(this));    % 6 x N
for k = 1:numel(this)
    xk = xhis(:,k);
    xdot(:,k) = f_fun(xk(1),xk(2),xk(3),xk(4),xk(5),xk(6));
end

Phi_tf = Phi(:,:,end);  % STM mapping epoch -> tf
disp('Phi(tf,t0) = ');
disp(Phi_tf);

% ======================  functions  ======================
function X_dot = EoMwithSTM(~, X, n, f, A)
  x   = X(1:n);                            
  Phi = reshape(X(n+1:n+n^2), n, n);       

  f1 = f(x(1),x(2),x(3),x(4),x(5),x(6));
  A  = A(x(1),x(2),x(3),x(4),x(5),x(6));

  dPhi = A*Phi;
  X_dot = [f1; dPhi(:)];

end

function [traj,Phi] = STM_compute(t0, tf, this, X0_aug, n, f_fun, A_fun, opts)
  traj = ode45(@(t,X) EoMwithSTM(t, X, n, f_fun, A_fun), [t0 tf], X0_aug, opts);

  Xhis = deval(traj, this);     
  N    = numel(this);

  Phi = zeros(n,n,N);
  for k = 1:N
      Phi(:,:,k) = reshape(Xhis(n+1:n+n^2, k), n, n);
  end
end
