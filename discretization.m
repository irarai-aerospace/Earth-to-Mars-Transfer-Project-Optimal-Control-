clc; clear;

dt = 1.0;
alpha = 0.5086;
g = [-3.7114e-3; 0; 0];

A = [zeros(3,3), eye(3),   zeros(3,1);
     zeros(3,3), zeros(3,3), zeros(3,1);
     zeros(1,3), zeros(1,3), 0];

B = [zeros(3,3), zeros(3,1);
     eye(3),     zeros(3,1);
     zeros(1,3), -alpha];

c = [zeros(3,1); g; 0];

Ak = expm(A*dt);

% Numerical integral for Bk, ck
nint = 2000;
tau = linspace(0, dt, nint);
Bk = zeros(7,4);
ck = zeros(7,1);

for i = 1:nint-1
    h = tau(i+1) - tau(i);
    tmid = 0.5*(tau(i) + tau(i+1));
    E = expm(A*(dt - tmid));
    Bk = Bk + E*B*h;
    ck = ck + E*c*h;
end

disp('Ak ='); disp(Ak);
disp('Bk ='); disp(Bk);
disp('ck ='); disp(ck);