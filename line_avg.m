%===============================================================================
% Calculates linear average between neighbouring elements of an array
%-------------------------------------------------------------------------------
function b = line_avg(d, w1, w2)

  if exist('w1','var')
    b = w1 .* d(1:end-1) + w2 .* d(2:end);
  else
    b = (d(2:end) + d(1:end-1))/2; 
  end

end
