%===============================================================================
% Discretizes momentum equation
% Units for system matrix for momentum are kg/s
%
%           |---o---|---o---|---o---|---o---|---o---|---o---|---o---|---o---|
%  cells:       1       2       3       4       5       6       7       8
%  faces:           1       2       3       4       5       6       7
%-------------------------------------------------------------------------------
function [a_u, t_u] = discretize_u(x_n, x_c, dx, dy, dz, dt, rho_c, mu_if);

  % Initialize sparse matrix
  a_u = diag(sparse(0));

  % Fetch the system size
  n_c = size(x_c, 2);

  %-------------------
  %
  % Inside the domain
  %
  %-------------------
  for c = 1:n_c

    %---------------
    % Unsteady term
    %---------------
    % Unit for a_u and t_u is: kg/m^3 * m^3 / s = kg/s
    t_u(c) = rho_c(c) * dx(c) * dy(c) * dz(c) / dt;
    a_u(c,c) = t_u(c);

    %---------------
    % Viscous terms
    %---------------
    % Unit for a_u is: kg/(m s) * m^2 / m = kg/s

    % This is the east side
    if(c < n_c)
      a = mu_if(c) * dy(c) * dz(c) / (x_c(c+1)-x_c(c));
      a_u(c,c+1) = -a;
      a_u(c,c)   =  a_u(c,c) + a;
    end

    % This is the west side
    if(c > 1)
      a = mu_if(c-1) * dy(c) * dz(c) / (x_c(c)-x_c(c-1));
      a_u(c,c-1) = -a;
      a_u(c,c)   =  a_u(c,c) + a;
    end

  end

  %-------------------
  %
  % On the boundaries
  %
  %-------------------
  a_u(1,1)     = a_u(1,1)     + dy(1)   * dz(1)   / (x_c(1)   - x_n(1));
  a_u(n_c,n_c) = a_u(n_c,n_c) + dy(n_c) * dz(n_c) / (x_c(n_c) - x_n(n_c));

end
