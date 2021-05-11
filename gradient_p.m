%===============================================================================
% This function calculates cell-centered gradients
%-------------------------------------------------------------------------------
%
%  Grid configuration:
%
%  nodes:   1       2       3       4       5       6       7       8       9
%           |---o---|---o---|---o---|---o---|---o---|---o---|---o---|---o---|
%  cells:       1       2       3       4       5       6       7       8
%  faces:           1       2       3       4       5       6       7
%
%-------------------------------------------------------------------------------
function g_p = gradient_p(x_n, x_c, w1, w2, p_c)

  % Fetch the size
  n_c = size(x_c, 2);

  % Inner face coordinates
  x_if = [x_n(2:n_c)];

  % Expand pressure to include boundary cells

  % Pressure at inner faces
  p_if = line_avg(p_c, w1, w2);  % size = [1, n_c-1]

  % Pressure at all faces (including boundary ones)
  p_af = [p_c(1), p_if, p_c(n_c)];  % size = [1, n_c+1]

  % Using excessive number of iterations (like in T-Flows)
  % in order to get a better comparison between two codes
  for i=1:40

    g_p = diff(p_af) ./ diff(x_n);

    for f = 1:n_c-1
      c1 = f;
      c2 = c1+1;
%     p_af(f+1) = w1(f) * p_c(c1) + w1(f) * (x_if(f) - x_c(c1)) * g_p(c1) ...
%               + w2(f) * p_c(c2) + w2(f) * (x_if(f) - x_c(c2)) * g_p(c2);
      p_af(f+1) = w1(f) * p_c(c1) + w2(f) * p_c(c2);
    end

    p_af(1)     = p_c(1)   - g_p(1)   * (x_c(1)     - x_n(1));
    p_af(n_c+1) = p_c(n_c) + g_p(n_c) * (x_n(n_c+1) - x_c(n_c));
  end

end
