%===============================================================================
% Performs simple linear interpolation of velocities at faces
%-------------------------------------------------------------------------------
function [u_if, u_af] = face_interpolation_01_rc_standard(  \
                        u_c, p_c, p_x,                      \
                        x_c, dv, a_u)

  % Unit for velocity m^3 * s / kg * kg/(m^2 s^2) = m / s

  u_if = line_avg(u_c + dv ./ spdiags(a_u, 0)' .* p_x);
  u_if = u_if - line_avg(dv ./ spdiags(a_u, 0)') .* diff(p_c) ./ diff(x_c);

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
