%===============================================================================
% Performs Rhie and Chow interpolation of momentum at faces, starting from
% interpolated velocities.  It is one in the series of functions:
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
function [u_if, u_af] = flux_03_rc_from_velocities(  ...
                        x_c, dv, a_u,                ...
                        u_c, p_c, p_x)

  % Form a helping array
  dv_au = dv ./ spdiags(a_u, 0)';

  %-------------------------------------------
  % Subtract cell-centered pressure gradients
  %-------------------------------------------

  % Unit for velocity: m^3 * s/kg * kg/(m^2 s^2) = m/s
  u_if = line_avg(u_c + dv_au .* p_x);                      % Rhie and Chow

   % Some info for comparison with T-Flows to print
   [ line_avg(u_c)',                   ...
     (diff(p_c))',                     ...
     (line_avg(dv_au))'                ...
     (line_avg(dv_au) .* diff(p_c))',  ...
     (line_avg(dv_au .* p_x))',        ...
     p_x(1:end-1)', p_x(2:end)' ]

  %----------------------------------
  % Add staggered pressure gradients
  %----------------------------------
  u_if = u_if - line_avg(dv_au) .* diff(p_c) ./ diff(x_c);  % Rhie and Chow

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
