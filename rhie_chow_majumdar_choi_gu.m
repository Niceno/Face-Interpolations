%===============================================================================
% Performs Rhie and Chow, Majumdar and Choi's interpolation of momentum
% at faces, starting from interpolated velocities.  It is one in the
% series of functions:
%
%   rhie_chow.m
%   rhie_chow_majumdar.m
%   rhie_chow_majumdar_choi.m
%   rhie_chow_majumdar_choi_gu.m
%   rhie_chow_choi.m
%   rhie_chow_choi_gu.m
%
% which all contain evolutionary steps Mencinger and Zun were doing to (7)
% and following equations in their paper in JCP from 2007.
%-------------------------------------------------------------------------------
function [u_if, u_af] = rhie_chow_majumdar_choi_gu(           ...
                        x_c, dv, urf_u, a_u, t_u, f_c, f_if,  ...
                        u_c, u_c_o, u_c_star,                 ...
                        u_if_o, u_if_star,                    ...
                        p_c, p_x)

  % Form helping arrays
  dv_au = dv  ./ spdiags(a_u, 0)';
  tu_au = t_u ./ spdiags(a_u, 0)';

  %------------------------------
  % Subtract cell-centered terms
  %------------------------------

  % Unit for velocity:
  % m^3 * s / kg * kg/(m^2 s^2) = m / s
  % m/s
  u_if = line_avg(u_c + dv_au  .* p_x               ...        % Rhie and Chow
                      - dv_au .* f_c                ...        % Choi
                      - tu_au .* u_c_o              ...        % Choi
                      - (1.0 - urf_u) * u_c_star);             % Majumdar

  %---------------------
  % Add staggered terms
  %---------------------
  u_if = u_if - line_avg(dv_au) .* diff(p_c) ./ diff(x_c) ...  % Rhie and Chow
              + line_avg(dv_au) .* f_if                   ...   % Gu
              + line_avg(tu_au) .* u_if_o                 ...  % Choi
              + (1.0 - urf_u) * u_if_star;                     % Majumdar

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
