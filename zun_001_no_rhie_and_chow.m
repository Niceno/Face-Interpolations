%-------------------------------------------------------------------------------
%
%  Grid configuration:
%
%  nodes:   1       2       3       4       5       6       7       8       9
%           |---o---|---o---|---o---|---o---|---o---|---o---|---o---|---o---|
%  cells:       1       2       3       4       5       6       7       8
%  faces:           1       2       3       4       5       6       7
%
%  Units:
%
%  velocity                [m/s]
%  pressure                [kg/(m s^2)]
%  gradients of pressure:  [kg/(m^2 s^2)]
%  density                 [kg/m^3]
%  dynamic viscosity:      [kg/(m s)]
%
%-------------------------------------------------------------------------------
clear;

g = 10.0;  % gravitational constant

%-------------------
% Solver parameters
%-------------------
par = load("-ascii", "solver.par")
n_steps = par(1);
dt      = par(2);
n_iters = par(3);
o_tol   = par(4);
u_iters = par(5);
urf_u   = par(6);
tol_u   = par(7);
p_iters = par(8);
urf_p   = par(9);
tol_p   = par(10);

%---------------------
% Physical properties
%---------------------
dat = load("-ascii", "physical.dat")
rho_water = dat(1);
rho_air   = dat(2);
mu_water  = dat(3);
mu_air    = dat(4);

%-----------------------------------------------------
% Read initial vof, it will also give number of cells
%-----------------------------------------------------
vof_c = load("-ascii", "vof.ini")';
n_c = size(vof_c, 2);

%-----------------
% Grid definition
%-----------------
l    =   0.1           % lenght of the domain
x_n  = l*(0:n_c)/n_c;  % node coordinates,       size = [1, n_c+1]
x_c  = line_avg(x_n);  % cell coordinates,       size = [1, n_c]
x_if = line_avg(x_c);  % inner face coordinates, size = [1, n_c-1]

dx = x_n(2:end) - x_n(1:end-1);  % size = [1, n_c]
dy = dx;                         % size = [1, n_c]
dz = dx;                         % size = [1, n_c]
dv = dx .* dy .* dz;             % size = [1, n_c]

%-------------------------------------------------------
% Distribution of physical properties on numerical mesh
%-------------------------------------------------------

% Values in cell centers
rho_c = vof_c * rho_air + (1-vof_c) * rho_water;  % density at cells
mu_c  = vof_c * mu_air  + (1-vof_c) * mu_water;   % visosity at cells

% Values in faces
% These are padded with boundary face values
rho_if = harm_avg(rho_c);                  % density at inner faces
mu_if  = harm_avg(mu_c);                   % viscosity at inner faces
rho_af = [rho_if(1),rho_if,rho_if(n_c-1)];   % append boundary values to all faces
mu_af  = [mu_if(1),mu_if,mu_if(n_c-1)];      % append boundary values to all faces

% Work out vof at faces 
% (If this is linear averaging, and rho_f is harmonic, 
% then the rho_f is disconnected from the vof at faces?)
vof_f = line_avg(vof_c);

% Work out the exact pressure solution
% (Interesting, but leave out for now)
% p_c(1:n_c) = 0;  % initialize pressure everywhere
% for i = 1:n_c-1
%   p_c(i+1) = p_c(i) + g * (x_c(i+1)-x_c(i)) * rho_f(i);
% end

%----------------
% Initial values
%----------------
u_c(1:n_c)  = 0.0;             % for velocities
p_c(1:n_c)  = 0.0;             % for pressure
pp_c(1:n_c) = 0.0;             % for pressure correction
p_x = gradient_p(x_c, p_c);  % for pressure gradients

%---------------------------------------
% Variables related to plotting results
%---------------------------------------
% fig_u = figure('Name', 'Velocity Iterations');
% fig_p = figure('Name', 'Pressure Iterations');
fig_a = figure('Name', 'Solution');

iter_leg = [];
step_leg = [];
iter_plot_int =   1;
step_plot_int = 100;
for i = 1:n_iters
  if(mod(i, iter_plot_int) == 0)
    iter_leg = [iter_leg; sprintf("Iter %d", i)];
  end
end
for i = 1:n_steps
  if(mod(i, step_plot_int) == 0)
    step_leg = [step_leg; sprintf("Step %d", i)];
  end
end

%----------------------
%
% Loop over time steps
%
%----------------------
for k = 1:n_steps

  % Store the last time step as old
  u_o = u_c;

  %----------------------------
  %
  % Loop over inner iterations
  %
  %----------------------------
  for iter = 1:n_iters

    %-------------------------------
    % Discretize momentum equations
    %-------------------------------
    % Unit for a_u is kg/s (see inside function for details)
    a_u = discretize_u(x_n, x_c, dx, dy, dz, rho_c, mu_if, dt);
    b_u(1:n_c) = 0.0;

    %------------------------------------------
    % Add unsteady term to the right hand side
    %------------------------------------------
    % Unit for unstead term: kg/(m^3) * m^3 * m / s / s = kg m/s^2
    b_u = b_u + rho_c .* dv .* u_o / dt;

    %------------------------------------------
    % Add pressure terms to momentum equations
    %------------------------------------------
    % Unit for presure term: kg/(m^2 s^2) m^3 = kg m/s^2
    b_u = b_u - p_x .* dv;

    %------------------------------------------------------
    % Cell centered buoyancy terms (that's wrong they say)
    %------------------------------------------------------
    % Unit for b_u is: m/s^2 * kg/m^3 * m^3 = kg m/s^2 = N
    b_u = b_u + g * rho_c .* dv;

    % Under-relax momentum equations
    for i=1:n_c
      b_u(i) = b_u(i) + (1.0-urf_u) / urf_u * a_u(i,i) * u_c(i);
      a_u(i,i) = a_u(i,i) / urf_u;
    end

    %--------------------------
    % Solve momentum equations
    %--------------------------
    % Units for velocity are: kg m/s^2 * s/kg = m/s
    u_c = pcg(a_u, b_u', tol_u, u_iters, [], [], u_c')';  % size = [1, n_c]

    if(exist('fig_u', 'var') == 1)
      plot_var(fig_u, 1, x_c, u_c, 'Cell Velocity Before Correction', iter);
    end

    %-------------------------------
    % Discretize pressure equations
    %-------------------------------
    % Units are: ms
    a_p = discretize_p(x_c, dy, dz, dv, a_u, rho_if);

    %--------------------------
    % Form r.h.s. for pressure
    %--------------------------
    % Unit for b_p is: kg/m^3 * m/s * m^2 = kg/s

    % Unit for velocity m/s
    u_if = line_avg(u_c);   % perform simple linear averaging (won't work)

    u_af = [0.0,u_if,0.0];  % append boundary values (just zeroes now)
    b_p = -diff(rho_af .* u_af) .* dy .* dz;

    if(exist('fig_u', 'var') == 1)
      plot_var(fig_u, 2, x_if, u_if, 'Face Velocity Before Correction', iter);
    end

    %-------------------------------
    % Solve for pressure correction
    %-------------------------------
    % Units for pressure are: kg/s / (ms) = kg/(m s^2)
    pp_c = pcg(a_p, b_p', tol_p, p_iters, [], [], pp_c')';  % size = [1, n_c]

    if(exist('fig_p', 'var') == 1)
      plot_var(fig_p, 1, x_c, pp_c, 'Pressure Correction', iter);
    end

    %---------------------------
    % Update the pressure field
    %---------------------------
    % Unit for pressure: N/m^2 = kg m/s^2 / m^2 = kg/(m s^2)
    p_c = p_c + pp_c * urf_p;

    if(exist('fig_p', 'var') == 1)
      plot_var(fig_p, 2, x_c, p_c, 'Pressure', iter);
    end

    %------------------------------------------------------
    % Calculate pressure and pressure correction gradients
    %------------------------------------------------------
    p_x  = gradient_p(x_c, p_c);
    pp_x = gradient_p(x_c, pp_c);

    if(exist('fig_p', 'var') == 1)
      plot_var(fig_p, 3, x_c, p_x,  'Pressure Gradient', iter);
      plot_var(fig_p, 4, x_c, pp_x, 'Pressure Correction Gradient', iter);
    end

    %-------------------------
    % Correct cell velocities
    %-------------------------
    u_c = u_c - pp_x .* dv ./ spdiags(a_u, 0)';

    if(exist('fig_u', 'var') == 1)
      plot_var(fig_u, 3, x_c, u_c, 'Cell Velocity After Correction', iter);
    end

    %-------------------------
    % Correct face velocities
    %-------------------------
    % Units for velocity: kg/(m s^2) * ms / m^2 * m^3 / kg = m/s
    for i=1:n_c-1
      a_f = 0.5 * (dy(i)*dz(i) + dy(i+1)*dz(i+1));
      u_if(i) = u_if(i) + (pp_c(i+1) - pp_c(i)) * a_p(i,i+1) / (rho_if(i) * a_f);
    end

    if(exist('fig_u', 'var') == 1)
      plot_var(fig_u, 4, x_if, u_if, 'Face Velocity After Correction', iter);
    end

  end

  o_res = sqrt(dot((u_c-u_o), (u_c-u_o)));
  printf("Outer loop residual %E\n", o_res);
  if(o_res < o_tol)
    ipoiopipoiop
  end

  if(exist('fig_a', 'var') == 1 && mod(k, step_plot_int) == 0)
    plot_var(fig_a, 1, x_c,  u_c,  'Cell Velocity',       k/step_plot_int);
    plot_var(fig_a, 2, x_if, u_if, 'Face Velocity',       k/step_plot_int);
    plot_var(fig_a, 3, x_c,  pp_c, 'Pressure Correction', k/step_plot_int);
    plot_var(fig_a, 4, x_c,  p_c,  'Pressure',            k/step_plot_int);
  end
end

% Place legends to plots
if(exist('fig_u', 'var') == 1)  figure(fig_u);  legend(iter_leg);  end
if(exist('fig_p', 'var') == 1)  figure(fig_p);  legend(iter_leg);  end

if(exist('fig_a', 'var') == 1)
  figure(fig_a);  subplot(2,2,1);  legend(step_leg);
end
