%===============================================================================
% Discretizes pressure equation from entries in momentum matrix
% Units for system matrix for pressure are ms
%
%           |---o---|---o---|---o---|---o---|---o---|---o---|---o---|---o---|
%  cells:       1       2       3       4       5       6       7       8
%  faces:           1       2       3       4       5       6       7
%-------------------------------------------------------------------------------
function a_p = discretize_p(x_c, dy, dz, dv, a_u, rho_if)

  % Initialize sparse matrix
  a_p = diag(sparse(0));

  % Fetch the system size
  n_c = size(a_u, 2);

  %-------------------
  % Inside the domain
  %-------------------
  for c = 1:n_c

    % Units: kg/m^3 * m * m^3 * s/kg = m/s

    % East side
    if(c < n_c)
      d = 0.5 * (dy(c)*dz(c) + dy(c+1)*dz(c+1)) / (x_c(c+1) - x_c(c));
      a = 0.5 * rho_if(c) * d * (dv(c) / a_u(c,c) + dv(c+1) / a_u(c+1, c+1));
      a_p(c,c+1) = -a;
      a_p(c,c)   =  a_p(c,c) + a;
    end

    % West side
    if(c > 1)
      d = 0.5 * (dy(c)*dz(c) + dy(c-1)*dz(c-1)) / (x_c(c) - x_c(c-1));
      a = 0.5 * rho_if(c-1) * d * (dv(c) / a_u(c,c) + dv(c-1) / a_u(c-1, c-1));
      a_p(c,c-1) = -a;
      a_p(c,c)   =  a_p(c,c) + a;
    endif

  end

end
