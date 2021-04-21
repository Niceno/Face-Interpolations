%===============================================================================
% This function calculates cell-centered gradients
%-------------------------------------------------------------------------------
function g_p = gradient_p(x_c, p_c)

  % Fetch the size
  n_c = size(x_c, 2);

  % First calculate pressure gradients at faces inside domain
  % These are correct, no question about it
  g_f = diff(p_c) ./ diff(x_c);   % size = [1, n_c-1]

  % Pad with boundary values by extrapolation
  % That's the best what you can do for a 2nd order method
  g_f = [g_f(1), g_f, g_f(n_c-1)];  % size = [1, n_c+1]

  % Finally take the average of face-gradients to get cell gradients
  % I hope this is analogue to least squares method in T-Flows
  g_p = line_avg(g_f);

end
