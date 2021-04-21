%===============================================================================
% Performs simple linear interpolation of velocities at faces
%-------------------------------------------------------------------------------
function [u_if, u_af] = face_interpolation_01_linear(u_c)

  % Unit for velocity m/s
  u_if = line_avg(u_c);   % perform simple linear averaging (won't work)

  u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)

end
