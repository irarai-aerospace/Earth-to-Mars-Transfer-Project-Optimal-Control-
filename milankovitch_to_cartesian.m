function [R_inertial, V_inertial] = milankovitch_to_cartesian(h,e,L,mu)
%h,e are column vector
  ecc = norm(e);
  hmag = norm(h);

  h_cap = h/norm(h);
  % e = e - (dot(e,h)/dot(h,h))*h;
  e_cap = e/norm(e);
  q_cap = cross(h_cap, e_cap);
  
  i = acos(h(3)/norm(h));
  
   % if abs(i) < 1e-12 || ecc < 1e-12
   %    error('Not defined');
   % end

  z_cap = [0;0;1];
  N = cross(z_cap,h);
  sig = atan2(N(2),N(1));
  om = atan2(dot(h_cap,cross(N,e)),dot(N,e));
  true_anomaly = mod(L - sig - om, 2*pi);
  p = hmag^2/mu;
  r = p/(1+ecc*cos(true_anomaly));
  R_inertial = r*((cos(true_anomaly))*e_cap + sin(true_anomaly)*q_cap);
  V_inertial = mu/hmag*(-sin(true_anomaly)*e_cap + (ecc+cos(true_anomaly))*q_cap);

end