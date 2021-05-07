%===============================================================================
% Conveniently plots a avariable from a simulation
%-------------------------------------------------------------------------------
function b = plot_cells(fig, sub, x, var, vof)

  % Switch to a desired figure ...
  figure(fig);

  % ... and the subplot
  subplot(2,2, sub);

  val = ylim;
  low_val = val(1)+ 0.333333 * (val(2) - val(1));
  hig_val = val(1)+ 0.666667 * (val(2) - val(1));

  % Plot lines representig cells
  for c=1:size(x,2)
    line( [x(c), x(c)], [val(1), val(2)], "Color", [200 200 200]/255,  ...
                                          "LineStyle", "-" );
  end

  % Plot bars for (future) vof
  for c=1:size(x,2)-1
    vof_val = low_val + vof(c) * (hig_val - low_val);

    line( [x(c), x(c+1)], [low_val, low_val], "Color", [200 200 200]/255,  ...
                                              "LineStyle", "-" );
    line( [x(c), x(c+1)], [hig_val, hig_val], "Color", [200 200 200]/255,  ...
                                              "LineStyle", "-" );
    line( [x(c), x(c+1)], [vof_val, vof_val], "Color", "c",      ...
                                              "LineStyle", "-",  ...
                                              "LineWidth", 2 );
  end

  for c=2:size(x,2)-1
    vof_val_m = low_val + vof(c-1) * (hig_val - low_val);
    vof_val_c = low_val + vof(c)   * (hig_val - low_val);
    line( [x(c), x(c)], [vof_val_m, vof_val_c], "Color", "c",    ...
                                              "LineStyle", "-",  ...
                                              "LineWidth", 2 );
  end

end
