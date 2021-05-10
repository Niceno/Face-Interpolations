%===============================================================================
% Performs Rhie and Chow, Choi's and Gu's interpolation of momentum at faces,
% starting from interpolated velocities.  It is one in the series of functions:
%
%   rhie_chow.m
%   rhie_chow_choi.m
%   rhie_chow_choi_gu.m
%
% which all contain evolutionary steps Mencinger and Zun were doing to (7)
% and following equations in their paper in JCP from 2007.
%-------------------------------------------------------------------------------
function [v_flux_if_n, v_flux_af_n] = rhie_chow_choi_gu(x_c, sx, dv,  ...
                                                        a_u, t_u,     ...
                                                        f_c, f_if,    ...
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

  %------------------------
  % Gu's correction part 1
  %------------------------
  u_c = u_c - v_m .* f_c;

  % Interpolate velocity
  u_f = line_avg(u_c);

  % Mimics bare-bone matrix from T-Flows
  % Unit: m
  a_fc = sx ./ diff(x_c);

  % Mimics pressure matrix from T-Flows
  % Unit: m^ 4 s / kg
  a12 = line_avg(v_m) .* a_fc;

  % Unit: 
  % m^3 s / kg * kg /(m^2 s^2) * m = m^2/s
  px_f = line_avg(v_m .* p_x) .* diff(x_c);

  %--------------------------
  % Rhie and Chow correction
  %--------------------------
  v_flux_if_n = u_f .* sx         ...
              - diff(p_c) .* a12  ...
              + px_f .* a_fc;

  %--------------------------
  % Choi's correction part 2
  %--------------------------
  v_flux_if_n = v_flux_if_n + line_avg(t_m) .* v_flux_o;

  %------------------------
  % Gu's correction part 2
  %------------------------
  v_flux_if_n = v_flux_if_n + line_avg(v_m) .* f_if .* sx;

  % Append boundary values (just zeroes now)
  v_flux_af_n = [0.0,v_flux_if_n,0.0];

end
