%===============================================================================
% Calculates linear average between neighbouring elements of an array
%-------------------------------------------------------------------------------
function b = weight_avg(d, w1, w2)

  b = w1 .* d(1:end-1) + w2 .* d(2:end);

end
