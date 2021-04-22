%===============================================================================
% Performs Rhie and Chow interpolation of velocities at faces
% This is a more elaborate version of flux_03_rc_standard_from_velocities,
% which follows closer the works on face iterpolation methods.
%
% It should contain evolutionary steps Mencinger and Zun were doing to (7)
% 
% 
%-------------------------------------------------------------------------------
function [u_if, u_af] = flux_08_rc_majumdar_choi_from_matrix(  ...
                        x_c, dv, urf_u, a_u, t_u, f_c,    ...
                        u_c, u_if_o, u_if_star, p_c, p_x)

  % Fetch the system size
  n_c = size(u_c, 2);

  % Calculate u_tilde (for equation (8) from Mencinger and Zun
  u_til(1:n_c) = 0.0;
  for c = 1:n_c
    if(c < n_c)  u_til(c) = u_til(c) - a_u(c,c+1) * u_c(c+1);  end  % east
    if(c > 1)    u_til(c) = u_til(c) - a_u(c,c-1) * u_c(c-1);  end  % west
    u_til(c) = u_til(c) / a_u(c,c);                                 % central
  end

  % Equation (10) from Mencinger and Zun
  % Unit for velocity:
  % m/s
  % kg/s / (kg/s) * m/s = m/s
  % m^3 / (kg/s) * N/m^3 = N s / kg = kg m/s^2 * s / kg = m/s
  % m/s
  % m^3 / (kg/s) * kg/(m^2 s^2) = m/s
  u_if = line_avg(  u_til                             ...  % (7.1) in my notes
                  + dv  ./ spdiags(a_u, 0)' .* f_c)   ...  % (7.4 and 7.5)
       + line_avg(t_u ./ spdiags(a_u, 0)') .* u_if_o  ...  % Choi
       + (1.0 - urf_u) * u_if_star                    ...  % Majumdar
       - line_avg(dv ./ spdiags(a_u, 0)') .* diff(p_c) ./ diff(x_c);

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
