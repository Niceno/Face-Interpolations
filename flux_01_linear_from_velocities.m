%===============================================================================
% Performs simple linear interpolation of velocities at faces
%
% This is as simple and as stupid as it gets, and it obviously, won't
% work with collocated variable arrangement
%-------------------------------------------------------------------------------
function [u_if, u_af] = flux_01_linear_from_velocities(u_c)

  % Unit for velocity m/s
  u_if = line_avg(u_c);   % perform simple linear averaging (won't work)

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
