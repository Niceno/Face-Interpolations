%===============================================================================
% Calculates linear average between neighbouring elements of an array
%-------------------------------------------------------------------------------
function b = line_avg(d)

  b = (d(2:end) + d(1:end-1))/2; 

end
