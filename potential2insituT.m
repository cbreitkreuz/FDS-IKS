function Temperature = potential2insituT(SALT, THETA, PRef, P)
% potential2insituT 
% Compute insitu temperature from potetnial temperatur.

% PRef = 0;
% P = abs(rC(k));

% C theta1
      del_P  = P - PRef;
      del_th = del_P .* sw_adtg(SALT,THETA,PRef);
      th     = THETA + 0.5 .* del_th;
      q      = del_th;
% C theta2
      del_th = del_P .* sw_adtg(SALT,th,PRef + 0.5 .* del_P);

      th     = th + (1 - 1/sqrt(2)).*(del_th - q);
      q      = (2-sqrt(2)).*del_th + (-2+3/sqrt(2)).*q;

% C theta3
      del_th = del_P .* sw_adtg(SALT,th,PRef + 0.5.*del_P);
      th     = th + (1 + 1/sqrt(2)) .* (del_th - q);
      q      = (2 + sqrt(2)) .* del_th + (-2-3/sqrt(2)) .* q;

% C theta4
      del_th = del_P.*sw_adtg(SALT,th,PRef+del_P);
      Temperature     = th + (del_th - 2 .* q)/(2*3);

      
end



%---------------------------------------------------------
function out = sw_adtg(SALT,THETA,P)

      sref = 35 .* 10^0;
      a0 =  3.5803 .*10^-5;
      a1 = 8.5258 .*10^-6;
      a2 = -6.836 .*10^-8;
      a3 =  6.6228 .*10^-10;

      b0 = 1.8932 .*10^-6;
      b1 = -4.2393 .*10^-8;

      c0 = 1.8741 .*10^-8;
      c1 = -6.7795 .*10^-10;
      c2 = 8.733 .*10^-12;
      c3 = -5.4481 .*10^-14;

      d0 = -1.1351 .*10^-10;
      d1 =  2.7759 .*10^-12;

      e0 = -4.6206 .*10^-13;
      e1 = 1.8676 .*10^-14;
      e2 = -2.1687 .*10^-16;

out =     a0 + (a1 + (a2 + a3.*THETA).*THETA).*THETA...
          + (b0 + b1 .* THETA).* (SALT-sref)...
          + ( (c0 + (c1 + (c2 + c3.*THETA).*THETA).*THETA) + (d0 + d1.*THETA).*(SALT-sref) ).*P...
          + (  e0 + (e1 + e2 .* THETA) .* THETA ) .* P .* P;
end
