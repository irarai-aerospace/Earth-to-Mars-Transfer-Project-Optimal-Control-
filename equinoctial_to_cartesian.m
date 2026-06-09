function [R_inertial, V_inertial] = equinoctial_to_cartesian(p,f,g,h,k,L,mu)
sig = atan2(k,h);
om = (atan2(g,f)) - sig;
i = 2*atan(sqrt(h^2+k^2));
nu = L - sig - om;
e = sqrt(f^2+g^2);
a = p/(1-e^2);
E  = 2*atan2( sqrt(1-e)*sin(nu/2), sqrt(1+e)*cos(nu/2) );
M = mod(E - e*sin(E), 2*pi);
[R_inertial, V_inertial] = keplerian_to_cartesian(a,e,M,i,sig,om,mu); %instantaneous M
end