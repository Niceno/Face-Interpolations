%===============================================================================
% Performs Rhie and Chow interpolation of momentum at faces, starting from
% interpolated matrix entries.  It is one in the series of functions:
%
%   flux_01_linear_from_velocities.m
%   flux_02_linear_from_matrix.m
%   flux_03_rc_from_velocities.m
%   flux_04_rc_from_matrix.m
%   flux_05_rc_majumdar_from_velocities.m
%   flux_06_rc_majumdar_from_matrix.m
%   flux_07_rc_majumdar_choi_from_velocities.m
%   flux_08_rc_majumdar_choi_from_matrix.m
%   flux_09_rc_majumdar_choi_gu_from_velocities.m
%   flux_10_rc_majumdar_choi_gu_from_matrix.m
%   flux_11_rc_choi_from_velocities.m
%   flux_12_rc_choi_from_matrix.m
%   flux_13_rc_choi_gu_from_velocities.m
%   flux_14_rc_choi_gu_from_matrix.m
%
% which all contain evolutionary steps Mencinger and Zun were doing to (7)
% and following equations in their paper in JCP from 2007.
%-------------------------------------------------------------------------------
function [u_if, u_af] = flux_04_rc_from_matrix(           ...
                        x_c, dv, urf_u, a_u, t_u, f_c,    ...
                        u_c, u_c_o, u_c_star, p_c, p_x)

  % Fetch the system size
  n_c = size(u_c, 2);

  % Calculate u_tilde (for equation (8) from Mencinger and Zun
  u_til(1:n_c) = 0.0;
  for c = 1:n_c
    if(c < n_c)  u_til(c) = u_til(c) - a_u(c,c+1) * u_c(c+1);  end  % east
    if(c > 1)    u_til(c) = u_til(c) - a_u(c,c-1) * u_c(c-1);  end  % west
    u_til(c) = u_til(c) / a_u(c,c);                                 % central
  end

  % Form helping arrays
  dv_au = dv  ./ spdiags(a_u, 0)';
  tu_au = t_u ./ spdiags(a_u, 0)';

  % Equation (10) from Mencinger and Zun
  % Unit for velocity:
  % m/s
  % kg/s / (kg/s) * m/s = m/s
  % m^3 / (kg/s) * N/m^3 = N s / kg = kg m/s^2 * s / kg = m/s
  % m/s
  % m^3 / (kg/s) * kg/(m^2 s^2) = m/s
  u_if = line_avg(  u_til                            ...   % (7.1) in my notes
                  + tu_au .* u_c_o                   ...   % (7.2 and 7.3)
                  + dv_au .* f_c                     ...   % (7.4 and 7.5)
                  + (1.0 - urf_u) * u_c_star)        ...   % (7.8)
       - line_avg(dv_au) .* diff(p_c) ./ diff(x_c);        % Rhie and Chow

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
