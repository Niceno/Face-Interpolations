%===============================================================================
% Performs linear interpolation of velocities at faces, but from matrix
% entries and righ hand side vectors.  This is an evolutionary step from
% flux_01_linear_from_velocities and it essentially serves to check if
% using the matrix entris (like in Peric's thesis and Mencinger's and Zun's
% paper from 2007) makes any difference in results.  Also, it is useful to
% check if all terms are properly discretized and used.
%
% This is equation (7) in paper from Mecinger and Zun (2007)
%-------------------------------------------------------------------------------
function [u_if, u_af] = flux_02_linear_from_matrix(  ...
                        dv, urf_u, a_u, t_u, f_c,    ...
                        u_c, u_c_o, p_x)

  % Fetch the system size
  n_c = size(u_c, 2);

  % Calculate u_tilde (for equation (8) from Mencinger and Zun
  u_til(1:n_c) = 0.0;
  for i = 1:n_c
    if(i < n_c)  u_til(i) = u_til(i) - a_u(i,i+1) * u_c(i+1);  end  % east
    if(i > 1)    u_til(i) = u_til(i) - a_u(i,i-1) * u_c(i-1);  end  % west
    u_til(i) = u_til(i) / a_u(i,i);                                 % central
  end

  % Equation (7) from Mencinger and Zun
  % Unit for velocity:
  % m/s
  % kg/s / (kg/s) * m/s = m/s
  % m^3 / (kg/s) * N/m^3 = N s / kg = kg m/s^2 * s / kg = m/s
  % m^3 / (kg/s) * kg/(m^2 s^2) = m/s
  u_if = line_avg(  u_til                             ...  % (7.1) in my notes
                  + t_u ./ spdiags(a_u, 0)' .* u_c_o  ...  % (7.2 and 7.3)
                  + dv  ./ spdiags(a_u, 0)' .* f_c    ...  % (7.4 and 7.5)
                  - dv  ./ spdiags(a_u, 0)' .* p_x    ...  % (7.6 and 7.7)
                  + (1.0 - urf_u) * u_c_o);                % (7.8)

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
