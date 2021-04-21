%===============================================================================
% Performs Rhie and Chow interpolation of velocities at faces
% This is how it was done in T-Flows when I started to work on surface tracking
%-------------------------------------------------------------------------------
function [u_if, u_af] = flux_03_rc_standard_from_velocities(  ...
                        x_c, dv, a_u,                         ...
                        u_c, p_c, p_x)

  % Unit for velocity m^3 * s / kg * kg/(m^2 s^2) = m / s

  u_if = line_avg(u_c + dv ./ spdiags(a_u, 0)' .* p_x);
  u_if = u_if - line_avg(dv ./ spdiags(a_u, 0)') .* diff(p_c) ./ diff(x_c);

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
