%===============================================================================
% Performs Rhie and Chow and Choi's interpolation of momentum at faces, 
% starting from interpolated velocities.  It is one in the series of functions:
%
%   rhie_chow.m
%   rhie_chow_choi.m
%   rhie_chow_choi_gu.m
%
% which all contain evolutionary steps Mencinger and Zun were doing to (7)
% and following equations in their paper in JCP from 2007.
%-------------------------------------------------------------------------------
function [v_flux_if_n, v_flux_af_n] = rhie_chow_choi(x_c, sx, dv,  ...
                                                     w1, w2,       ...
                                                     a_u, t_u,     ...
                                                     u_n, u_o,     ...
                                                     v_flux_o,     ...
                                                     p_c, p_x)

  % Take velocities as computed
  u_c = u_n;

  % Form helping arrays
  % Unit: m^3 s / kg
  v_m = dv  ./ spdiags(a_u, 0)';

  % Helping array for Choi
  t_m = t_u ./ spdiags(a_u, 0)';

  %--------------------------
  % Choi's correction part 1
  %--------------------------
  u_c = u_c - t_m .* u_o;

  % Interpolate velocity
  u_f = weight_avg(u_c, w1, w2);

  % Mimics bare-bone matrix from T-Flows
  % Unit: m
  a_fc = sx ./ diff(x_c);

  % Mimics pressure matrix from T-Flows
  % Unit: m^ 4 s / kg
  a12 = weight_avg(v_m, w1, w2) .* a_fc;

  % Unit: 
  % m^3 s / kg * kg /(m^2 s^2) * m = m^2/s
  px_f = weight_avg(v_m .* p_x, w1, w2) .* diff(x_c);

  %--------------------------
  % Rhie and Chow correction
  %--------------------------
  v_flux_if_n = u_f .* sx         ...
              - diff(p_c) .* a12  ...
              + px_f .* a_fc;

  %--------------------------
  % Choi's correction part 2
  %--------------------------
  v_flux_if_n = v_flux_if_n + weight_avg(t_m, w1, w2) .* v_flux_o;

  % Append boundary values (just zeroes now)
  v_flux_af_n = [0.0,v_flux_if_n,0.0];

end
