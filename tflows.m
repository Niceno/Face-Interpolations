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
%     u_n, u_o, v_flux_if_n, v_flux_af_n, p_c
%-------------------------------------------------------------------------------
clear

g        =  0.0;     % gravitational constant  [m/s^2]
mass_src =  3.0;     % mass source             [kg/s]
kappa    =  0.0;     % curvature               [1/m]
d_vof    =  0.0333;  % for moving interface

u_west = 0.0;
u_east = 0.0;

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
fid   = fopen('algorithms.def','r');
algor = fgetl(fid)
fclose(fid);

%---------------------
% Physical properties
%---------------------
dat = load('-ascii', 'physical.dat')
rho_water = dat(1);
rho_steam  = dat(2);
mu_water   = dat(3);
mu_steam   = dat(4);
sigma      = dat(5);

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

% Stretch the grid
d_x(1:n_c) = 1;
stretch = 1.0;
for i = 2:floor(n_c/2)
  d_x(i) = d_x(i-1) * stretch;
  d_x(n_c+1-i) = d_x(i);
end
d_x = d_x / sum(d_x) * l;
for i = 1:n_c
  x_n(i+1) = x_n(i) + d_x(i);
end

% Compute cell and face coordinates
x_c           = line_avg(x_n);  % cell coordinates,       size = [1, n_c]
x_if(1:n_c-1) = x_n(2:n_c);     % inner face coordinates, size = [1, n_c-1]

% Face-weigths; size = [1, n_c-1]
w1 = (x_c(2:n_c) - x_n(2:n_c)) ./ diff(x_c);
w2 = 1.0 - w1;

% test gradients: p_c = x_c;
% test gradients: p_x = gradient_p(x_n, x_c, w1, w2, p_c);
% test gradients: rtyuio

dx = x_n(2:end) - x_n(1:end-1);  % size = [1, n_c]
dy = dx * 2000.0;                % size = [1, n_c]
dz = dx * 2000.0;                % size = [1, n_c]
dv = dx .* dy .* dz;             % size = [1, n_c]
sx = dy(1) * dz(1)               % size = [1, n_c-1]

%----------------
% Initial values
%----------------
u_n(1:n_c)  = 0.0;                         % for velocities
v_flux_if_n = line_avg(u_n, w1, w2) * sx;  % vol. flux at inner faces [m^3/s]
p_c(1:n_c)  = 0.0;                         % for pressure
pp_c(1:n_c) = 0.0;                         % for pressure correction
p_x = gradient_p(x_n, x_c, w1, w2, p_c);   % for pressure gradients

%---------------------------------------
% Variables related to plotting results
%---------------------------------------

% Figures for iterations
iter_leg = [];
if(iter_plot_int > 0)
  fig_u = figure('Name', 'Velocity Iterations');
  fig_p = figure('Name', 'Pressure Iterations');
end
if(step_plot_int > 0)
  fig_a = figure('Name', ['Transients With ', strrep(algor, '_', ' '),  ...;
                          ' and dt=', mat2str(dt)]);
end

%----------------------
%
% Loop over time steps
%
%----------------------

% Initialize array for residuals in time step
o_res = [];

for step = 1:n_steps

  printf('#======================\n');
  printf('#                      \n');
  printf('# Time Step: %d        \n', step);
  printf('#                      \n');
  printf('#======================\n');

   %-------------
   % Advance VOF
   %-------------
   for f = 1:n_c-1
     c1 = f;
     c2 = c1+1;
     if(vof_c(c1) <= 0.5  && vof_c(c2) > 0)
       vof_c(c2) = vof_c(c2) - d_vof;
       vof_c(c2) = max(0, vof_c(c2));
       break;
     end
   end

  %------------------------------------------------------------------------
  % Distribution of physical properties on numerical mesh depending on vof
  %------------------------------------------------------------------------

  % Set some initial values in the cell centers ("_c").
  % For density, this will be used for inertial term.
  % (See the next comment)
  rho_c = vof_c * rho_water + (1-vof_c) * rho_steam;  % density at cells
  mu_c  = vof_c * mu_water  + (1-vof_c) * mu_steam;   % visosity at cells

  % Work out the values in the faces.  Here you can use linear or harmonic mean
  % You could have, in fact, started from prescribing physical properties from
  % faces, which would make sense only if you initially prescribed vof in faces.
  % (I am not sure, maybe it is something worth considering, could be.)
  rho_if = harm_avg(rho_c);                    % density at inner faces

  mu_if  = harm_avg(mu_c);                     % viscosity at inner faces
  mu_af  = [mu_if(1),mu_if,mu_if(n_c-1)];      % append boundary values

  % Density-biased height coefficients
  % (The application so far wasn't successful.  It helps to smooth out
  %  velocities for Stefan's case, but ruins all the others.)
  h1 = 0.5 * rho_if ./ rho_c(1:n_c-1);
  h2 = 0.5 * rho_if ./ rho_c(2:n_c);

  %----------------------------------------------
  % Store the last time step as old (suffix "o")
  %----------------------------------------------
  u_o  = u_n;
  v_flux_o = v_flux_if_n;

  %----------------------------
  %
  % Loop over inner iterations
  %
  %----------------------------
  for iter = 1:n_iters

    %-----------------------------------------
    % Initialize right-hand side for momentum
    %-----------------------------------------
    % Unit for a_u is kg/s (see inside function for details)
    b_u(1:n_c) = 0.0;

    %-------------------------------
    % Discretize momentum equations
    %-------------------------------
    [a_u, t_u, b_u] = discretize_u(x_n, x_c, dx, sx, dt,  ...
                                   rho_c, mu_af,          ...
                                   u_west, u_east);

    % Store original matrix like you would do in T-Flows
    m_u = a_u';

    %--------------------------------------
    % Discretize pressure equations before
    %  momentum system was under-relaxed
    %--------------------------------------
    % Units are: m^4 s / kg
    a_p = discretize_p(x_c, sx, dv, w1, w2, a_u);

    % Under-relax the discretized momentum equations
    for c=1:n_c
      a_u(c,c) = a_u(c,c) / urf_u;
    end

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

    %-------------------------------------
    % Face-centered surface tension terms
    %-------------------------------------
    f_if = -sigma * kappa * diff(vof_c) ./ diff(x_c);  % N/m * 1/m * 1/m = N/m^3

    %------------------------------
    % Face centered buoyancy terms
    %------------------------------
    % Unit for b_u is: m/s^2 * kg/m^3 * m^3 = kg m/s^2 = N
    f_if = f_if + g * rho_if;   % kg/m^3 * m/s^2 = kg/(m^2 s^2) = N/m^3

    % Face force due to mass transfer
    % kg^2/s^2 * m^3/kg / m^2 / m^3 = kg / (m^2 s^2) = N / m^3
    % These lead to nothing, they were way
    % too small compared to pressure forces
    f_mt_if(1:n_c-1) = 0.0;
    for f = 1:n_c-1
      c1 = f;
      c2 = c1+1;
      if(vof_c(c1) <= 0.5  && vof_c(c2) > 0.5)
        f_mt_if(f) = mass_src^2 * (1/rho_water+1/rho_steam) ...
                   / sx / (0.5*(dv(c1)+dv(c2)));
      end
    end
    f_mt_af = [f_mt_if(1), f_mt_if, f_mt_if(n_c-1)];   % N
    f_mt_c = line_avg(f_mt_af);

    % Expand forces to boundary faces too
    f_af = [f_if(1), f_if, f_if(n_c-1)];               % N/m^3
    % made no difference :-( f_af = f_af + f_mt_af;

    %-----------------------------------------------------------------------
    %
    %   This line is absolutelly crucial.  It works out density at cells
    %   (because they will determing also the forces in cells) from facial
    %   values.  This interpolation must be linear - I think because of
    %   the linear assumptions made in the finite volume method.
    %
    %   The suffix "_f" means this density is used to compute facial forces
    %
    %-----------------------------------------------------------------------
    f_c = line_avg(f_af);   % kg/m^3 * m/s^2 = kg/(m^2 s^2) = N/m^3
    b_u = b_u + f_c .* dv;

    % Under-relax forces in momentum equations
    % (a_u is already divided by urf_u above)
    for c=1:n_c
      b_u(c) = b_u(c) + (1.0-urf_u) * a_u(c,c) * u_n(c);
    end

    % Store initial residual for this time step
    % (I use the same name as Mencinger and Zun)
    if(iter == 1) r_n_0 = norm(a_u * u_n' - b_u'); end
    if(step == 1) r_1_0 = r_n_0; end

    %--------------------------
    % Solve momentum equations
    %--------------------------

    % Units for velocity are: kg m/s^2 * s/kg = m/s
    u_n = pcg(a_u, b_u', tol_u, u_iters, [], [], u_n')';  % size = [1, n_c]

    if(mod(iter, iter_plot_int) == 0)
      iter_leg = [iter_leg; sprintf('iter %d', iter)];
      plot_var(fig_u, 1, x_c, u_n, 'Cell Velocity Before Correction',  ...
                                   iter / iter_plot_int);
    end

    %--------------------------
    % Form r.h.s. for pressure
    %--------------------------

    % Perform interpolation at a cell face
    switch(algor)
      case 'Rhie-Chow'
        [v_flux_if_n, v_flux_af_n] = rhie_chow(x_c, sx, dv, ...
                                               w1, w2,      ...
                                               m_u,         ...
                                               u_n,         ...
                                               p_c, p_x);
      case 'Rhie-Chow_Choi'
        [v_flux_if_n, v_flux_af_n] = rhie_chow_choi(x_c, sx, dv,  ...
                                                    w1, w2,       ...
                                                    m_u, t_u,     ...
                                                    u_n, u_o,     ...
                                                    v_flux_o,     ...
                                                    p_c, p_x);
      case 'Rhie-Chow_Choi_Gu'
        [v_flux_if_n, v_flux_af_n] = rhie_chow_choi_gu(x_c, sx, dv,  ...
                                                       w1, w2,       ...
                                                       m_u, t_u,     ...
                                                       f_c, f_if,    ...
                                                       u_n, u_o,     ...
                                                       v_flux_o,     ...
                                                       p_c, p_x);
      otherwise
        do_something_completely_different ();
    end

    % Unit for b_p is: m^3/s
    b_p = -diff(v_flux_af_n);

    % Add volumetric source (boiling)
    for f = 1:n_c-1
      c1 = f;
      c2 = c1+1;
      if(vof_c(c1) <= 0.5  && vof_c(c2) > 0.5)
        vol_src_steam =  mass_src / rho_steam;  % volume of steam created
        vol_src_water = -mass_src / rho_water;  % volume of water lost
        b_p(c1)  = b_p(c1) + vol_src_water;
        b_p(c2)  = b_p(c2) + vol_src_steam;
        b_p(n_c) = b_p(n_c) - (vol_src_steam + vol_src_water);

        % Mimic convective outflow
        u_east = mass_src / sx;
      end
    end

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_u, 2, x_if, v_flux_if_n, 'Face Flux Before Correction',  ...
                                            iter / iter_plot_int);
    end

    %-------------------------------
    % Solve for pressure correction
    %-------------------------------
    % Units for pressure are: kg/s / (ms) = kg/(m s^2)
    pp_c = pcg(a_p, b_p', tol_p, p_iters, [], [], pp_c')';  % size = [1, n_c]

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_p, 1, x_c, pp_c, 'Pressure Correction',  ...
                                    iter / iter_plot_int);
    end

    %---------------------------
    % Update the pressure field
    %---------------------------
    % Unit for pressure: N/m^2 = kg m/s^2 / m^2 = kg/(m s^2)
    p_c = p_c + pp_c * urf_p;

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_p, 2, x_c, p_c, 'Pressure',  ...
                                   iter / iter_plot_int);
    end

    %------------------------------------------------------
    % Calculate pressure and pressure correction gradients
    %------------------------------------------------------
    p_x  = gradient_p(x_n, x_c, w1, w2, p_c);
    pp_x = gradient_p(x_n, x_c, w1, w2, pp_c);

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_p, 3, x_c, pp_x, 'Pressure Correction Gradient',  ...
                                   iter / iter_plot_int);
      plot_var(fig_p, 4, x_c, p_x,  'Pressure Gradient',  ...
                                    iter / iter_plot_int);
    end

    %-------------------------
    % Correct cell velocities
    %-------------------------
    u_n = u_n - pp_x .* dv ./ spdiags(m_u, 0)';

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_u, 3, x_c, u_n, 'Cell Velocity After Correction',  ...
                                   iter / iter_plot_int);
    end

    %-------------------------
    % Correct face velocities
    %-------------------------
    % Units for velocity: kg/(m s^2) * ms / m^2 * m^3 / kg = m/s
    for c=1:n_c-1
      v_flux_if_n(c) = v_flux_if_n(c) + (pp_c(c+1) - pp_c(c)) * a_p(c,c+1);
    end

    if(mod(iter, iter_plot_int) == 0)
      plot_var(fig_u, 4, x_if, v_flux_if_n, 'Face Flux After Correction',  ...
                                            iter/iter_plot_int);
    end

    % Work out residuals in the current iteration
    % (I use the same name as Mencinger and Zun)
    r_n_i = norm(a_u * u_n' - b_u')

    res_dt = r_n_i / r_n_0;
    printf('Residual reduction: %E\n', res_dt);

    if(res_dt < eps_dt)
      printf('Convergence in time step reached\n');
      break;
    end

  end


  % Plot transient solutions
  if(mod(step, step_plot_int) == 0)
    plot_var(fig_a, 1, x_c,  u_n,         'Cell Velocity',       step/step_plot_int);
    plot_var(fig_a, 2, x_if, v_flux_if_n, 'Face Flux',           step/step_plot_int);
    plot_var(fig_a, 3, x_c,  pp_c,        'Pressure Correction', step/step_plot_int);
    plot_var(fig_a, 4, x_c,  p_c,         'Pressure',            step/step_plot_int);
  end

  % Check if steady state has been reached
  res_st = r_n_0 / r_1_0;
  printf('Residual to steady state: %E\n', res_st);
  o_res = [o_res, res_st];
  if(res_st < eps_st)
    printf('#=====================================\n');
    printf('#                                     \n');
    printf('# Steady State Reached In %d Steps \n', step);
    printf('#                                     \n');
    printf('#=====================================\n');
    break;
  end

end

% Plot cells in all graphs for transient solution
if(step_plot_int > 0)
  plot_cells(fig_a, 1, x_n,  u_n,         vof_c);
  plot_cells(fig_a, 2, x_n,  v_flux_if_n, vof_c);
  plot_cells(fig_a, 3, x_n,  pp_c,        vof_c);
  plot_cells(fig_a, 4, x_n,  p_c,         vof_c);

  step_leg = [];
  for c = 1:step
    if(mod(c, step_plot_int) == 0)
      step_leg = [step_leg; sprintf('step %d', c)];
    end
  end
  figure(fig_a);  subplot(2,2,2);  legend(step_leg);
end

% Place legends to plots for iterations
if(iter_plot_int > 0)
  plot_cells(fig_u, 1, x_n, u_n,         vof_c);
  plot_cells(fig_u, 3, x_n, u_n,         vof_c);
  plot_cells(fig_u, 2, x_n, v_flux_if_n, vof_c);
  plot_cells(fig_u, 4, x_n, v_flux_if_n, vof_c);

  plot_cells(fig_p, 1, x_n, pp_c, vof_c);
  plot_cells(fig_p, 2, x_n, p_c,  vof_c);
  plot_cells(fig_p, 3, x_n, pp_x, vof_c);
  plot_cells(fig_p, 4, x_n, p_x,  vof_c);

  figure(fig_u);  legend(iter_leg);
  figure(fig_p);  legend(iter_leg);
end

% Plot final solution
fig_f = figure('Name', ['Final Solution With ', strrep(algor, '_', ' '), ...
                        ' and dt=', mat2str(dt)]);
plot_var(fig_f, 1, x_c,  u_n,         'Cell Velocity',       1);
plot_var(fig_f, 2, x_if, v_flux_if_n, 'Face Flux',           1);
plot_var(fig_f, 3, x_c,  pp_c,        'Pressure Correction', 1);
plot_var(fig_f, 4, x_c,  p_c,         'Pressure',            1);
plot_cells(fig_f, 1, x_n,  u_n,         vof_c);
plot_cells(fig_f, 2, x_n,  v_flux_if_n, vof_c);
plot_cells(fig_f, 3, x_n,  pp_c,        vof_c);
plot_cells(fig_f, 4, x_n,  p_c,         vof_c);

% Plot residual history
if(n_steps > 1)
  figure('Name', 'Residual History');
  x=[]; for c=1:size(o_res,2) x=[x, c]; end
  semilogy(x, o_res, '-*');
  title(['Convergence History With ', strrep(algor, '_', ' '), ...
         ' and dt=', mat2str(dt)]);
end
