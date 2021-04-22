%===============================================================================
% Performs Rhie and Chow interpolation of velocities at faces
% This is how it was done in T-Flows when I started to work on surface tracking
%-------------------------------------------------------------------------------
function [u_if, u_af] = flux_07_rc_majumdar_choi_from_velocities(  ...
                        x_c, dv, urf_u, a_u, t_u,                  ...
                        u_c, u_c_o, u_c_star,                      ...
                        u_if_o, u_if_star,                         ...
                        p_c, p_x)

  %------------------------------
  % Subtract cell-centered terms
  %------------------------------

  % Unit for velocity:
  % m^3 * s / kg * kg/(m^2 s^2) = m / s
  % m/s
  u_if = line_avg(u_c + dv ./ spdiags(a_u, 0)'  .* p_x    ...  % Rhie and Chow
                      - t_u ./ spdiags(a_u, 0)' .* u_c_o  ...  % Choi
                      - (1.0 - urf_u) * u_c_star);             % Majumdar

  %---------------------
  % Add staggered terms
  %---------------------
  u_if = u_if - line_avg(dv ./ spdiags(a_u, 0)') .* diff(p_c) ./ diff(x_c) ...
       + line_avg(t_u ./ spdiags(a_u, 0)') .* u_if_o  ...  % Choi
       + (1.0 - urf_u) * u_if_star;                        % Majumdar

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
