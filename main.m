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
%
%  Calling convention:
%
%  When calling other scripts, it would be good to have some conventions on
%  the order of arguments, for the sake of consistency.  I might order them
%  in terms of their actual accuracy. I might do like this:
%
%  1. Geometrical parameters (as firmly defined) from the essential to the
%     derived ones, hence:
%
%     x_n, x_c, x_if, dx, dy, dz, dv
%
%  2. Time step, also very clearly defined
%
%     dt
%
%  3. Physical properties (which are less certain than physical dimensions,
%     and are somehow interpolated on the grid), in the order:
%
%     rho_c, rho_if, rho_af, mu_c, mu_if, mu_af
%
%  4. Under-relaxation factors
%
%     urf_u, urf_p
%
%  5. Matrices and right hand sides (since they depend on physical properties
%     and discretization methods) in the order:
%
%     a_u, t_u, a_p, f_c, b_u, b_p
%
%  6. Unknowns, as the least accurate, come last
%
%     u_c, u_c_o, u_if, u_af, p_c
%-------------------------------------------------------------------------------
clear;

g = 10.0;  % gravitational constant

%-------------------
% Solver parameters
%-------------------
par = load("-ascii", "solver.par")
dt            = par(1);   % time step
n_steps       = par(2);   % number of time steps
step_plot_int = par(3);   % plot interval for time steps
n_iters       = par(4);   % number of simple iterations
iter_plot_int = par(5);   % plot interval for iterations
o_tol         = par(6);   % outer (time step) loop tolerance
urf_u         = par(7);   % under-relaxation factor for momentum
u_iters       = par(8);   % number of iterations for momentum solver
tol_u         = par(9);   % tolerance for momentum solver
urf_p         = par(10);  % under-relaxation factor for pressure
p_iters       = par(11);  % number of iterations for pressure solver
tol_p         = par(12);  % tolerance for pressure solver

%------------------------------
% Face interpolation algorithm
%------------------------------
fid = fopen('algorithms.def','r');
algor = fgetl(fid);
fclose(fid);

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

% Figures for iterations
iter_leg = [];
if(iter_plot_int > 0)
  fig_u = figure('Name', 'Velocity Iterations');
  fig_p = figure('Name', 'Pressure Iterations');
  for i = 1:n_iters
    if(mod(i, iter_plot_int) == 0)
      iter_leg = [iter_leg; sprintf("iter %d", i)];
    end
  end
end

% Figures for time steps
step_leg = [];
if(step_plot_int > 0)
  fig_a = figure('Name', ['Transients With ', algor]);
  fig_f = figure('Name', ['Final Solution With ', algor]);
  for i = 1:n_steps
    if(mod(i, step_plot_int) == 0)
      step_leg = [step_leg; sprintf("step %d", i)];
    end
  end
end

% Initialize array for residuals in time step
o_res = [];

%----------------------
%
% Loop over time steps
%
%----------------------
for k = 1:n_steps

  printf("#========================\n");
  printf("#                        \n");
  printf("# Time Step: %d          \n", k);
  printf("#                        \n");
  printf("#========================\n");
  % Store the last time step as old

  u_c_o = u_c;

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
    [a_u, t_u] = discretize_u(x_n, x_c, dx, dy, dz, dt, rho_c, mu_if);
    b_u(1:n_c) = 0.0;

    %------------------------------------------
    % Add unsteady term to the right hand side
    %------------------------------------------
    % Unit for unstead term: kg/(m^3) * m^3 * m / s / s = kg m/s^2
    b_u = b_u + rho_c .* dv .* u_c_o / dt;

    %------------------------------------------
    % Add pressure terms to momentum equations
    %------------------------------------------
    % Unit for presure term: kg/(m^2 s^2) m^3 = kg m/s^2
    b_u = b_u - p_x .* dv;

    %------------------------------------------------------
    % Cell centered buoyancy terms (that's wrong they say)
    %------------------------------------------------------
    % Unit for b_u is: m/s^2 * kg/m^3 * m^3 = kg m/s^2 = N
    f_c = g * rho_c;        % kg/m^3 * m/s^2 = kg/(m^2 s^2) = N/m^3
    b_u = b_u + f_c .* dv;

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

    if(exist('fig_u', 'var') == 1 && mod(iter, iter_plot_int) == 0)
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

    % Perform interpolation at a cell face
    switch(algor)
      case 'Linear_From_Velocities'
        [u_if, u_af] = flux_01_linear_from_velocities(u_c);
      case 'Linear_From_Matrix'
        [u_if, u_af] = flux_02_linear_from_matrix(           ...
                       dv, urf_u, a_u, t_u, f_c,             ...
                       u_c, u_c_o, p_x);
      case 'RC_Standard_From_Velocities'
        [u_if, u_af] = flux_03_rc_standard_from_velocities(  ...
                       x_c, dv, a_u,                         ...
                       u_c, p_c, p_x);
      case 'RC_Standard_From_Matrix'
        [u_if, u_af] = flux_04_rc_standard_from_matrix(      ...
                       x_c, dv, urf_u, a_u, t_u, f_c,        ...
                       u_c, u_c_o, p_c, p_x);
      otherwise
        do_something_completely_different ();
    end

    % Unit for b_p is: kg/m^3 * m/s * m^2 = kg/s
    b_p = -diff(rho_af .* u_af) .* dy .* dz;

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_u, 2, x_if, u_if, 'Face Velocity Before Correction', iter);
    end

    %-------------------------------
    % Solve for pressure correction
    %-------------------------------
    % Units for pressure are: kg/s / (ms) = kg/(m s^2)
    pp_c = pcg(a_p, b_p', tol_p, p_iters, [], [], pp_c')';  % size = [1, n_c]

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_p, 1, x_c, pp_c, 'Pressure Correction', iter);
    end

    %---------------------------
    % Update the pressure field
    %---------------------------
    % Unit for pressure: N/m^2 = kg m/s^2 / m^2 = kg/(m s^2)
    p_c = p_c + pp_c * urf_p;

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_p, 2, x_c, p_c, 'Pressure', iter);
    end

    %------------------------------------------------------
    % Calculate pressure and pressure correction gradients
    %------------------------------------------------------
    p_x  = gradient_p(x_c, p_c);
    pp_x = gradient_p(x_c, pp_c);

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_p, 3, x_c, p_x,  'Pressure Gradient', iter);
      plot_var(fig_p, 4, x_c, pp_x, 'Pressure Correction Gradient', iter);
    end

    %-------------------------
    % Correct cell velocities
    %-------------------------
    u_c = u_c - pp_x .* dv ./ spdiags(a_u, 0)';

    if(mod(iter, iter_plot_int) == 0)
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

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_u, 4, x_if, u_if, 'Face Velocity After Correction', iter);
    end

  end

  o_res = [o_res, sqrt(dot((u_c-u_c_o), (u_c-u_c_o)))];
  printf("Outer loop residual %E\n", o_res(end));
  if(o_res(end) < o_tol)
    break;
  end

  % Plot transient solutions
  if(mod(k, step_plot_int) == 0)
    plot_var(fig_a, 1, x_c,  u_c,  'Cell Velocity',       k/step_plot_int);
    plot_var(fig_a, 2, x_if, u_if, 'Face Velocity',       k/step_plot_int);
    plot_var(fig_a, 3, x_c,  pp_c, 'Pressure Correction', k/step_plot_int);
    plot_var(fig_a, 4, x_c,  p_c,  'Pressure',            k/step_plot_int);
  end
end

% Plot final solution
plot_var(fig_f, 1, x_c,  u_c,  'Cell Velocity',       1);
plot_var(fig_f, 2, x_if, u_if, 'Face Velocity',       1);
plot_var(fig_f, 3, x_c,  pp_c, 'Pressure Correction', 1);
plot_var(fig_f, 4, x_c,  p_c,  'Pressure',            1);

% Place legends to plots
if(mod(iter, iter_plot_int) == 0)
  figure(fig_u);  legend(iter_leg);
  figure(fig_p);  legend(iter_leg);
end

if(step_plot_int > 0)
  figure(fig_a);  subplot(2,2,1);  legend(step_leg);
end
