%===============================================================================
% Discretizes pressure equation from entries in momentum matrix
% Units for system matrix for pressure are ms
%
%           |---o---|---o---|---o---|---o---|---o---|---o---|---o---|---o---|
%  cells:       1       2       3       4       5       6       7       8
%  faces:           1       2       3       4       5       6       7
%-------------------------------------------------------------------------------
function a_p = discretize_p(x_c, sx, dv, w1, w2, a_u)

  % Initialize sparse matrix
  a_p = diag(sparse(0));

  % Fetch the system size
  n_c = size(a_u, 2);

  %-------------------
  % Inside the domain
  %-------------------
  for f = 1:n_c-1

    c1 = f;
    c2 = c1 + 1;

    % Units: m * m^3 * s/kg = m^4 s / kg
    d = sx / (x_c(c2) - x_c(c1));
    a = w1(f) * d * (dv(c1) / a_u(c1, c1)) ...
      + w2(f) * d * (dv(c2) / a_u(c2, c2));
    a_p(c1,c2) = -a;
    a_p(c2,c1) = -a;
    a_p(c1,c1) =  a_p(c1,c1) + a;
    a_p(c2,c2) =  a_p(c2,c2) + a;

  end

end
