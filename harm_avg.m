%===============================================================================
% Calculates harmonic average between neighbouring elements of an array
%-------------------------------------------------------------------------------
function b = harm_avg(d)

    b = 2.0 ./ (1.0 ./ d(2:end) + 1.0 ./ d(1:end-1));

end
