%===============================================================================
% Performs Rhie and Chow interpolation of momentum at faces, starting from
% interpolated velocities.  It is one in the series of functions:
%
%   rhie_chow.m
%   rhie_chow_choi.m
%   rhie_chow_choi_gu.m
%
% which all contain evolutionary steps Mencinger and Zun were doing to (7)
% and following equations in their paper in JCP from 2007.
%-------------------------------------------------------------------------------
function [v_flux_if_n, v_flux_af_n] = rhie_chow(x_c, sx, dv,  ...
                                                w1, w2,       ...
                                                a_u,          ...
                                                u_n,          ...
                                                p_c, p_x)

  % Take velocities as computed
  u_c = u_n;

  % Form helping arrays
  % Unit: m^3 s / kg
  v_m = dv  ./ spdiags(a_u, 0)';

  % Interpolate velocity
  u_f = line_avg(u_c, w1, w2);

  % Mimics bare-bone matrix from T-Flows
  % Unit: m
  a_fc = sx ./ diff(x_c);

  % Mimics pressure matrix from T-Flows
  % Unit: m^ 4 s / kg
  a12 = line_avg(v_m, w1, w2) .* a_fc;

  % Unit: 
  % m^3 s / kg * kg /(m^2 s^2) * m = m^2/s
  px_f = line_avg(v_m .* p_x, w1, w2) .* diff(x_c);

  %--------------------------
  % Rhie and Chow correction
  %--------------------------
  v_flux_if_n = u_f .* sx         ...
              - diff(p_c) .* a12  ...
              + px_f .* a_fc;

  % Append boundary values (just zeroes now)
  v_flux_af_n = [0.0,v_flux_if_n,0.0];

end
