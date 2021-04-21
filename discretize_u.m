%===============================================================================
% Discretizes momentum equation
% Units for system matrix for momentum are kg/s
%
%           |---o---|---o---|---o---|---o---|---o---|---o---|---o---|---o---|
%  cells:       1       2       3       4       5       6       7       8
%  faces:           1       2       3       4       5       6       7
%-------------------------------------------------------------------------------
function a_u = discretize_u(x_n, x_c, dx, dy, dz, rho_c, mu_if, dt)

  % Initialize sparse matrix
  a_u = diag(sparse(0));

  % Fetch the system size
  n_c = size(x_c, 2);

  % Inside the domain
  for i = 1:n_c

    % Unsteady term
    % Unit for a_u is: kg/m^3 * m^3 / s = kg/s
    a_u(i,i) = rho_c(i) * dx(i) * dy(i) * dz(i) / dt;

    % Viscous terms
    % Unit for a_u is: kg/(m s) * m^2 / m = kg/s
    if(i < n_c)
      a = mu_if(i) * dy(i) * dz(i) / (x_c(i+1)-x_c(i));
      a_u(i,i+1) = -a;
      a_u(i,i)   =  a_u(i,i) + a;
    end

    if(i > 1)
      a = mu_if(i-1) * dy(i) * dz(i) / (x_c(i)-x_c(i-1));
      a_u(i,i-1) = -a;
      a_u(i,i)   =  a_u(i,i) + a;
    end

  end

  % On the boundaries
  a_u(1,1)     = a_u(1,1)     + dy(1)   * dz(1)   / (x_c(1)   - x_n(1));
  a_u(n_c,n_c) = a_u(n_c,n_c) + dy(n_c) * dz(n_c) / (x_c(n_c) - x_n(n_c));

end
