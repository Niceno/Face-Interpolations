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
  for i = 1:n_c

    % Units: kg/m^3 * m * m^3 * s/kg = m/s

    % East side
    if(i < n_c)
      d = 0.5 * (dy(i)*dz(i) + dy(i+1)*dz(i+1)) / (x_c(i+1) - x_c(i));
      a = 0.5 * rho_if(i) * d * (dv(i) / a_u(i,i) + dv(i+1) / a_u(i+1, i+1));
      a_p(i,i+1) = -a;
      a_p(i,i)   =  a_p(i,i) + a;
    end

    % West side
    if(i > 1)
      d = 0.5 * (dy(i)*dz(i) + dy(i-1)*dz(i-1)) / (x_c(i) - x_c(i-1));
      a = 0.5 * rho_if(i-1) * d * (dv(i) / a_u(i,i) + dv(i-1) / a_u(i-1, i-1));
      a_p(i,i-1) = -a;
      a_p(i,i)   =  a_p(i,i) + a;
    endif

  end

end
