%===============================================================================
% This function calculates cell-centered gradients
%-------------------------------------------------------------------------------
function g_p = gradient_p(x_n, x_c, p_c)

  % Fetch the size
  n_c = size(x_c, 2);

  % Expand coordinates to include boundary cells
  x_e = [x_n(1), x_c, x_n(n_c+1)];

  % Expand pressure to include boundary cells
  p_e = [p_c(1), p_c, p_c(n_c)];  % size = [1, n_c+2]

  % Using excessive number of iterations (like in T-Flows)
  % in order to get a better comparison between two codes
  for i=1:36

    % First calculate pressure gradients at faces inside and on the boundary
    % These are correct, no question about it
    g_f = diff(p_e) ./ diff(x_e);   % size = [1, n_c+1]

    % Finaly, get the cell centered average
    g_p = line_avg(g_f);            % size = [1, n_c]

    p_e(1)     = p_c(1)   - g_p(1)   * (x_c(1)     - x_n(1));
    p_e(n_c+2) = p_c(n_c) + g_p(n_c) * (x_n(n_c+1) - x_c(n_c));

  end

end
