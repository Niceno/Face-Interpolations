%===============================================================================
% Conveniently plots a avariable from a simulation
%-------------------------------------------------------------------------------
function b = plot_var(fig, sub, x, p, name, iter)

  % Switch to a desired figure ...
  figure(fig);

  % ... and the subplot
  subplot(2,2, sub);

  % Plot what you want ...
  plot(x, p, '-*');

  % ... and give it a title
  title(name);

  % Hold the plot if you are here for the first time
  if(iter == 1)
    hold;
  end

end
