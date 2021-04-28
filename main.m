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
%     a_u, t_u, a_p, f_c, f_if, b_u, b_p
%
%  6. Unknowns, as the least accurate, come last
%
%     u_c, u_c_o, u_c_star, u_if, u_af, p_c
%-------------------------------------------------------------------------------

g = 10.0;  % gravitational constant

%-------------------
% Solver parameters
%-------------------
par = load('-ascii', 'solver.par')
dt            = par( 1);  % time step
n_steps       = par( 2);  % number of time steps
step_plot_int = par( 3);  % plot interval for time steps
eps_st        = par( 4);  % steady state tolerance
n_iters       = par( 5);  % number of simple iterations
iter_plot_int = par( 6);  % plot interval for iterations
eps_dt        = par( 7);  % outer (time step) loop tolerance
urf_u         = par( 8);  % under-relaxation factor for momentum
u_iters       = par( 9);  % number of iterations for momentum solver
tol_u         = par(10);  % tolerance for momentum solver
urf_p         = par(11);  % under-relaxation factor for pressure
p_iters       = par(12);  % number of iterations for pressure solver
tol_p         = par(13);  % tolerance for pressure solver

%------------------------------
% Face interpolation algorithm
%------------------------------
fid = fopen('algorithms.def','r');
algor   = fgetl(fid)
p_matrix= fgetl(fid)
fclose(fid);

%---------------------
% Physical properties
%---------------------
dat = load('-ascii', 'physical.dat')
rho_water = dat(1);
rho_air   = dat(2);
mu_water  = dat(3);
mu_air    = dat(4);

%-----------------------------------------------------
% Read initial vof, it will also give number of cells
%-----------------------------------------------------
vof_c = load('-ascii', 'vof.ini')';
n_c = size(vof_c, 2);

%-----------------
% Grid definition
%-----------------
l    =   0.1           % lenght of the domain
x_n  = l*(0:n_c)/n_c;  % node coordinates,       size = [1, n_c+1]
x_c  = line_avg(x_n);  % cell coordinates,       size = [1, n_c]
x_if = line_avg(x_c);  % inner face coordinates, size = [1, n_c-1]

dx = x_n(2:end) - x_n(1:end-1);  % size = [1, n_c]
dy = dx * 2000.0;                % size = [1, n_c]
dz = dx * 2000.0;                % size = [1, n_c]
dv = dx .* dy .* dz;             % size = [1, n_c]

%-------------------------------------------------------
% Distribution of physical properties on numerical mesh
%-------------------------------------------------------

% Set some "initial" values in the cell centers
% (See the next comment)
rho_c = vof_c * rho_air + (1-vof_c) * rho_water;  % density at cells
mu_c  = vof_c * mu_air  + (1-vof_c) * mu_water;   % visosity at cells

% Work out the values in the faces.  Here you can use linear or harmonic mean
% You could have, in fact, started from prescribing physical properties from
% faces, which would make sense only if you initially prescribed vof in faces.
% (I am not sure, maybe it is something worth considering, could be.)
rho_if = harm_avg(rho_c);                    % density at inner faces
rho_af = [rho_if(1),rho_if,rho_if(n_c-1)];   % append boundary values

mu_if  = harm_avg(mu_c);                     % viscosity at inner faces
mu_af  = [mu_if(1),mu_if,mu_if(n_c-1)];      % append boundary values

% This line is absolutelly crucial.  It works out density at cells (because
% they will determing also the forces in cells) from facial values.  This
% interpolation must be linear - I think because of the linear assumptions
% made in the finite volume method.
% rho_c  = line_avg(rho_af);                   % recompute from face values

% Work out vof at faces 
% (If this is linear averaging, and rho_f is harmonic, 
% then the rho_f is disconnected from the vof at faces?)
vof_f = line_avg(vof_c);

% Work out the exact pressure solution
% (Interesting, but leave out for now)
% p_c(1:n_c) = 0;  % initialize pressure everywhere
% for c = 1:n_c-1
%   p_c(c+1) = p_c(c) + g * (x_c(c+1)-x_c(c)) * rho_f(c);
% end

%----------------
% Initial values
%----------------
u_c(1:n_c)  = 0.0;             % for velocities
u_if = line_avg(u_c);
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
  for c = 1:n_iters
    if(mod(c, iter_plot_int) == 0)
      iter_leg = [iter_leg; sprintf('iter %d', c)];
    end
  end
end
if(step_plot_int > 0)
  fig_a = figure('Name', ['Transients With ', strrep(algor, '_', ' '),  ...;
                          ' and dt=', mat2str(dt)]);
end

%----------------------
%
% Discretize equations
%
%----------------------

%-------------------------------
% Discretize momentum equations
%-------------------------------
[a_u, t_u] = discretize_u(x_n, x_c, dx, dy, dz, dt, rho_c, mu_af);

%---------------------------------------------------------
% Discretize pressure equations if it is not over-relaxed
% meaning it is formed before momentum is under-relaxed
%----------------------------------------------------------
if(strcmp(p_matrix, 'Pressure_Matrix_Default'))
  % Units are: ms
  a_p = discretize_p(x_c, dy, dz, dv, a_u);
end

% Under-relax the discretized momentum equations
for c=1:n_c
  a_u(c,c) = a_u(c,c) / urf_u;
end

%------------------------------------------------------
% Discretize pressure equations if it is over-relaxed
% meaning it is formed after momentum is under-relaxed
%------------------------------------------------------
if(! strcmp(p_matrix, 'Pressure_Matrix_Default'))
  % Units are: ms
  a_p = discretize_p(x_c, dy, dz, dv, a_u);
end

%-------------------------------
% Discretize pressure equations
%-------------------------------

%----------------------
%
% Loop over time steps
%
%----------------------

% Initialize array for residuals in time step
o_res = [];

for k = 1:n_steps

  printf('#======================\n');
  printf('#                      \n');
  printf('# Time Step: %d        \n', k);
  printf('#                      \n');
  printf('#======================\n');

  % Store the last time step as old (suffix "o")
  u_c_o  = u_c;
  u_if_o = u_if;

  %----------------------------
  %
  % Loop over inner iterations
  %
  %----------------------------
  for i = 1:n_iters

    % Store velocity from last iteration (suffix "star")
    u_c_star  = u_c;
    u_if_star = u_if;

    %-----------------------------------------
    % Initialize right-hand side for momentum
    %-----------------------------------------
    % Unit for a_u is kg/s (see inside function for details)
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
    f_if = g * rho_if;   % kg/m^3 * m/s^2 = kg/(m^2 s^2) = N/m^3
    f_c  = g * rho_c;    % kg/m^3 * m/s^2 = kg/(m^2 s^2) = N/m^3
    b_u  = b_u + f_c .* dv;

    % Under-relax forces in momentum equations
    % (a_u is already divided by urf_u above)
    for c=1:n_c
      b_u(c) = b_u(c) + (1.0-urf_u) * a_u(c,c) * u_c(c);
    end

    % Store initial residual for this time step
    % (I use the same name as Mencinger and Zun)
    if(i == 1) r_n_0 = norm(a_u * u_c' - b_u'); end
    if(k == 1) r_1_0 = r_n_0; end

    %--------------------------
    % Solve momentum equations
    %--------------------------

    % Store the previous iteration (suffix "p" for previous)
    u_c_p = u_c;

    % Units for velocity are: kg m/s^2 * s/kg = m/s
    u_c = pcg(a_u, b_u', tol_u, u_iters, [], [], u_c')';  % size = [1, n_c]

    u_c_before_correction = u_c;

    if(exist('fig_u', 'var') == 1 && mod(i, iter_plot_int) == 0)
      plot_var(fig_u, 1, x_c, u_c, 'Cell Velocity Before Correction', i);
    end

    %--------------------------
    % Form r.h.s. for pressure
    %--------------------------

    % Perform interpolation at a cell face
    switch(algor)
      case 'Linear_From_Velocities'
        [u_if, u_af] = flux_01_linear_from_velocities(u_c);
      case 'Linear_From_Matrix'
        [u_if, u_af] = flux_02_linear_from_matrix(                   ...
                       dv, urf_u, a_u, t_u, f_c,                     ...
                       u_c, u_c_o, u_c_star, p_x);
      case 'Rhie-Chow_Standard_From_Velocities'
        [u_if, u_af] = flux_03_rc_standard_from_velocities(          ...
                       x_c, dv, a_u,                                 ...
                       u_c, p_c, p_x);
      case 'Rhie-Chow_Standard_From_Matrix'
        [u_if, u_af] = flux_04_rc_standard_from_matrix(              ...
                       x_c, dv, urf_u, a_u, t_u, f_c,                ...
                       u_c, u_c_o, u_c_star, p_c, p_x);
      case 'Rhie-Chow_Majumdar_From_Velocities'
        [u_if, u_af] = flux_05_rc_majumdar_from_velocities(          ...
                       x_c, dv, urf_u, a_u,                          ...
                       u_c, u_c_star, u_if_star, p_c, p_x);
      case 'Rhie-Chow_Majumdar_From_Matrix'
        [u_if, u_af] = flux_06_rc_majumdar_from_matrix(              ...
                       x_c, dv, urf_u, a_u, t_u, f_c,                ...
                       u_c, u_c_o, u_if_star, p_c, p_x);
      case 'Rhie-Chow_Majumdar_Choi_From_Velocities'
        [u_if, u_af] = flux_07_rc_majumdar_choi_from_velocities(     ...
                       x_c, dv, urf_u, a_u, t_u,                     ...
                       u_c, u_c_o, u_c_star,                         ...
                       u_if_o, u_if_star,                            ...
                       p_c, p_x);
      case 'Rhie-Chow_Majumdar_Choi_From_Matrix'
        [u_if, u_af] = flux_08_rc_majumdar_choi_from_matrix(         ...
                       x_c, dv, urf_u, a_u, t_u, f_c,                ...
                       u_c, u_if_o, u_if_star, p_c, p_x);
      case 'Rhie-Chow_Majumdar_Choi_Gu_From_Velocities'
        [u_if, u_af] = flux_09_rc_majumdar_choi_gu_from_velocities(  ...
                       x_c, dv, urf_u, a_u, t_u, f_c, f_if,          ...
                       u_c, u_c_o, u_c_star,                         ...
                       u_if_o, u_if_star,                            ...
                       p_c, p_x);
      case 'Rhie-Chow_Majumdar_Choi_Gu_From_Matrix'
        [u_if, u_af] = flux_10_rc_majumdar_choi_gu_from_matrix(      ...
                       x_c, dv, urf_u, a_u, t_u, f_if,               ...
                       u_c, u_if_o, u_if_star, p_c, p_x);
      otherwise
        do_something_completely_different ();
    end

    % Unit for b_p is: m/s * m^2 = m^3/s
    b_p = -diff(u_af) .* dy .* dz;

    if(mod(i, iter_plot_int) == 0)
      plot_var(fig_u, 2, x_if, u_if, 'Face Velocity Before Correction', i);
    end

    %-------------------------------
    % Solve for pressure correction
    %-------------------------------
    % Units for pressure are: kg/s / (ms) = kg/(m s^2)
    pp_c = pcg(a_p, b_p', tol_p, p_iters, [], [], pp_c')';  % size = [1, n_c]

    if(mod(i, iter_plot_int) == 0)
      plot_var(fig_p, 1, x_c, pp_c, 'Pressure Correction', i);
    end

    %---------------------------
    % Update the pressure field
    %---------------------------
    % Unit for pressure: N/m^2 = kg m/s^2 / m^2 = kg/(m s^2)
    p_c = p_c + pp_c * urf_p;

    if(mod(i, iter_plot_int) == 0)
      plot_var(fig_p, 2, x_c, p_c, 'Pressure', i);
    end

    %------------------------------------------------------
    % Calculate pressure and pressure correction gradients
    %------------------------------------------------------
    p_x  = gradient_p(x_c, p_c);
    pp_x = gradient_p(x_c, pp_c);

    if(mod(i, iter_plot_int) == 0)
      plot_var(fig_p, 3, x_c, p_x,  'Pressure Gradient', i);
      plot_var(fig_p, 4, x_c, pp_x, 'Pressure Correction Gradient', i);
    end

    %-------------------------
    % Correct cell velocities
    %-------------------------
    u_c = u_c - pp_x .* dv ./ spdiags(a_u, 0)';

    if(mod(i, iter_plot_int) == 0)
      plot_var(fig_u, 3, x_c, u_c, 'Cell Velocity After Correction', i);
    end

    %-------------------------
    % Correct face velocities
    %-------------------------
    % Units for velocity: kg/(m s^2) * ms / m^2 * m^3 / kg = m/s
    for c=1:n_c-1
      a_f = 0.5 * (dy(c)*dz(c) + dy(c+1)*dz(c+1));
      u_if(c) = u_if(c) + (pp_c(c+1) - pp_c(c)) * a_p(c,c+1) / (rho_if(c) * a_f);
    end

    if(mod(i, iter_plot_int) == 0)
      plot_var(fig_u, 4, x_if, u_if, 'Face Velocity After Correction', i);
    end

    % Work out residuals in the current iteration
    % (I use the same name as Mencinger and Zun)
    r_n_i = norm(a_u * u_c' - b_u')

    res_dt = r_n_i / r_n_0;
    printf('Residual reduction: %E\n', res_dt);

    if(res_dt < eps_dt)
      printf('Convergence in time step reached\n');
      break;
    end

  end


  % Plot transient solutions
  if(mod(k, step_plot_int) == 0)
    plot_var(fig_a, 1, x_c,  u_c,  'Cell Velocity',       k/step_plot_int);
    plot_var(fig_a, 2, x_if, u_if, 'Face Velocity',       k/step_plot_int);
    plot_var(fig_a, 3, x_c,  pp_c, 'Pressure Correction', k/step_plot_int);
    plot_var(fig_a, 4, x_c,  p_c,  'Pressure',            k/step_plot_int);
  end

  % Check if steady state has been reached
  res_st = r_n_0 / r_1_0;
  printf('Residual to steady state: %E\n', res_st);
  o_res = [o_res, res_st];
  if(res_st < eps_st)
    printf('#======================\n');
    printf('#                      \n');
    printf('# Steady State Reached \n');
    printf('#                      \n');
    printf('#======================\n');
    break;
  end

end

% Plot final solution
fig_f = figure('Name', ['Final Solution With ', strrep(algor, '_', ' '), ...
                        ' and dt=', mat2str(dt)]);
plot_var(fig_f, 1, x_c,  u_c,  'Cell Velocity',       1);
plot_var(fig_f, 2, x_if, u_if, 'Face Velocity',       1);
plot_var(fig_f, 3, x_c,  pp_c, 'Pressure Correction', 1);
plot_var(fig_f, 4, x_c,  p_c,  'Pressure',            1);

% Plot residual history
figure('Name', 'Residual History');
x=[]; for c=1:size(o_res,2) x=[x, c]; end
semilogy(x, o_res, '-*');
title(['Convergence History With ', strrep(algor, '_', ' '), ...
       ' and dt=', mat2str(dt)]);

% Place legends to plots
if(mod(i, iter_plot_int) == 0)
  figure(fig_u);  legend(iter_leg);
  figure(fig_p);  legend(iter_leg);
end

if(step_plot_int > 0)
  step_leg = [];
  for c = 1:k
    if(mod(c, step_plot_int) == 0)
      step_leg = [step_leg; sprintf('step %d', c)];
    end
  end
  figure(fig_a);  subplot(2,2,1);  legend(step_leg);
end

