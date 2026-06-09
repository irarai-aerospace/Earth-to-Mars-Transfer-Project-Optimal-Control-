function [R_inertial, V_inertial] = keplerian_to_cartesian(a,e,M,i,sig,om,mu)

  [E,~] = Kepler_equation(e,M);

  true_anomaly = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2)) ;

   rc = a*(1-e*cos(E));
  
  A = [cos(true_anomaly);sin(true_anomaly);0];
  r = rc*A;

  B = [-sin(E);sqrt(1-e^2)*cos(E);0];
  K = sqrt(mu*a)/rc;
  v = K*B;
  Rotation = rotz(rad2deg(sig)) * rotx(rad2deg(i)) * rotz(rad2deg(om));
  
  R_inertial = Rotation*r;
  V_inertial = Rotation*v;

end



