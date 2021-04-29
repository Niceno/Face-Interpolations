%===============================================================================
% Performs linear interpolation of momentum at faces, starting from 
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
function [u_if, u_af] = flux_01_linear_from_velocities(u_c)

  % Unit for velocity m/s
  u_if = line_avg(u_c);   % perform simple linear averaging (won't work)

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
