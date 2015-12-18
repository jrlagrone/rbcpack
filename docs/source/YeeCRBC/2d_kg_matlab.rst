.. highlight:: matlab

*****************************
2D Klein-Gordon with a Corner
*****************************

This tutorial explains the code `klein_gordon_2d_corner_DAB.m <https://bitbucket.org/rbcpack/rbcpack/src/default/YeeCRBC/Matlab/klein_gordon_2d_corner_DAB.m>`_ and requires the use of `optimal_cosines.m <https://bitbucket.org/rbcpack/rbcpack/src/default/YeeCRBC/Matlab/optimal_cosines.m>`_

Introduction
============

This Matlab code implements a second order finite difference approximation to
the 2D Klein-Gordon equation. On two intersecting sides, the grid is terminated 
with a Double Absorbing Boundary (DAB).

The primary thing to notice here is that the DAB is essentially identical to the
1D case described in the :doc:`1D Klein-Gordon example <1d_kg_matlab>` on each
of the individual sides. At the corner where the two DAB sides intersect, the 
procedure is again similar, but now we have a double indexed set of auxiliary
variables.

What this program does
----------------------

For this example, we consider the 2D Klein-Gordon equation,

.. math::
  :nowrap:

  \begin{align}
    \frac{d^2 u}{d t^2} &= c^2 \left(\frac{d^2 u}{d x^2} + \frac{d^2 u}{d y^2} \right) - su,
  \end{align}

where :math:`c>0`.

For this example, we will impose Dirichlet boundary conditions on the left side
in the x-direction and at the bottom in the y-direction. We truncate the domain 
at the top and the right side with a DAB. We will use zero initial conditions.

.. math::
  :nowrap:

  \begin{align}
    u(x,y,0) = 0, \\
    \frac{du}{dt}(x,y,0) = 0, \\
    u(0,y,t) = 0, \\
    u(x,0,t) = 0.
  \end{align}

.. _discretization:

We discretize this equation using second order finite differences on the 
domain :math:`[0, .1]` with mesh spacing of :math:`\Delta x = \Delta y` and we define

.. math::
  :nowrap:

  \begin{align}
    x_i &=  i\Delta x, \qquad i=0,...,n-1,
    y_i &=  i\Delta y, \qquad i=0,...,n-1,
  \end{align}

We choose a time step size, :math:`\Delta t`, satisfying

.. math::
  :nowrap:

  \begin{align}
    \Delta t \leq \frac{dx}{c \sqrt{2}}.
  \end{align}

Letting 

.. math::
  :nowrap:

  \begin{align}
    t_n = n \Delta t,
  \end{align}

:math:`u` is on the grid by:

.. math::
  :nowrap:

  \begin{align}
    \frac{u(x_i,y_j,t_{n-1}) - 2u(x_i,y_j,t_{n}) + u(x_i,y_j,t_{n+1})}{\Delta t} 
    = c\left( \frac{u(x_{i-1},y_j,t_{n}) - 2u(x_i,y_j,t_{n}) + u(x_{i+1},y_j,t_{n})}{\Delta x} 
      + \frac{u(x_{i},y_{j-1},t_{n}) - 2u(x_i,y_j,t_{n}) + u(x_{i},y_{j+1},t_{n})}{\Delta y}\right)
      - s u(x_i,y_j,t_{n}).
  \end{align}

We use the discretization of the DAB described in :ref:`theory overview <dab_discretization>`


We drive the simulation with a point source that takes the form of
a differentiated Gaussian.


The commented program
=====================

We begin by choosing some basic simulation parameters. First we choose the number
of grid points to use in the discetization. Then we choose the number of grid
points to extend the domain by so we can compare to a larger simulation to check
error. In this case, we use 250 grid points and extend the larger simulation by
250 grid points, which corresponds to running on the domain [0,0.2]. Then, we
choose the problem and source parameters. ::

  n = 250;           % number of grid points in domain [0,0.1]
  m = 250;           % number of grid points to extent the domain by for a reference
                     % solution using a larger simulation
  c = 1;             % wave speed
  s = 20;            % dispersion parameter
  nsteps = 1000;     % number of time steps
  cfl = 0.99; % cfl ratio (1 is exact dispersion relation, but num. unstable)

  % plotting
  figure(1)
  plot_freq = 10;   % plot frequency (in time steps)

  % compute grid spacing / time step
  h = 0.1 / (n - 1); % we'll just use dx=dy
  dt = cfl / (c*sqrt(2/h^2));

  % source paramters
  tw = 25*dt;             % pulse width
  t0 = 5*tw;              % pulse delay (roughly need t0 > 4*tw to be smooth)
  amp = 1;                % pulse "amplitude"
  sloc = [60, 25];      % source location

Next, we choose the DAB parameters. We choose the number of recursions to use,
:math:`p` and how wide the DAB layer should be in grid points. We require at 
least three points to support the update to the Klein-Gordon equation. For
efficiency and accuracy, 3 is the best choice; however, if the auxiliary variables
are to be plotted, increasing the thickness is desireable.

In this example, we opt to use the optimal cosines. ::

  % choose optimal cosines
  T = nsteps * dt;        % simulation length
  delta = max(n-sloc)*h; % Seperation of DAB from source
  eta = delta / (c * T);
  if (p>0)
    [at, errest] = optimal_cosines(eta, p-1);
    a = at(1:2:2*p);
    ab = at(2:2:2*p);
    sig = 0.5*dt*(1 - a.*a) ./ (T*a);
    sigb = 0.5*dt*(1 - ab.*ab) ./ (T*ab);
  end

  fprintf('The reflection coefficient is %e \n', errest);

Next, we allocate the storage for all of the field values and auxiliary variables
we will use. ::

  % allocate storage
  unew = zeros(n); % field values
  ucur = zeros(n);
  uold = zeros(n);

  runew = zeros(n+m); % for larger reference simulation
  rucur = zeros(n+m);
  ruold = zeros(n+m);

  udabnew_east = zeros(ndab, n, p+1); % east dab aux. variables
  udabcur_east = zeros(ndab, n, p+1);
  udabold_east = zeros(ndab, n, p+1);

  udabnew_north = zeros(n, ndab, p+1); % north dab aux. variables
  udabcur_north = zeros(n, ndab, p+1);
  udabold_north = zeros(n, ndab, p+1);

  udabnew_ne = zeros(ndab, ndab, p+1, p+1); % northeast corner dab aux. variables
  udabcur_ne = zeros(ndab, ndab, p+1, p+1);
  udabold_ne = zeros(ndab, ndab, p+1, p+1);

We begin time stepping by updating all of the internal field values and adding a
source term. ::

  % time step
  for t=1:nsteps
    
    % internal updates --- eqn 54, in DAB paper
    unew(2:n-1, 2:n-1) = 2*ucur(2:n-1, 2:n-1) - uold(2:n-1, 2:n-1) ...
        + ((c*dt)/h)^2 * (ucur(3:n, 2:n-1) - 4*ucur(2:n-1,2:n-1) ...
        + ucur(1:n-2, 2:n-1) + ucur(2:n-1, 1:n-2) + ucur(2:n-1, 3:n)) ...
        - s*dt^2*ucur(2:n-1, 2:n-1);
    
    % reference solution (same thing on a bigger domain)
    runew(2:n+m-1, 2:n+m-1) = 2*rucur(2:n+m-1, 2:n+m-1) - ruold(2:n+m-1, 2:n+m-1) ...
        + ((c*dt)/h)^2 * (rucur(3:n+m, 2:n+m-1) - 4*rucur(2:n+m-1,2:n+m-1) ...
        + rucur(1:n+m-2, 2:n+m-1) + rucur(2:n+m-1, 1:n+m-2) + rucur(2:n+m-1, 3:n+m)) ...
        - s*dt^2*rucur(2:n+m-1, 2:n+m-1);
    
    % add source
    unew(sloc(1), sloc(2)) = unew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    runew(sloc(1), sloc(2)) = runew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);

Next we update the DAB layer at the top of the domain. To begin the DAB update, 
in all auxiliary variables we use the same update equation that we use to evolve
the interior points. ::
    
    %% North DAB
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    udabnew_north(2:n-1, 2:ndab-1,:) = 2*udabcur_north(2:n-1, 2:ndab-1,:) - udabold_north(2:n-1, 2:ndab-1,:) ...
        + ((c*dt)/h)^2 * (udabcur_north(3:n, 2:ndab-1,:) - 4*udabcur_north(2:n-1,2:ndab-1,:) ...
        + udabcur_north(1:n-2, 2:ndab-1,:) + udabcur_north(2:n-1, 1:ndab-2,:) + udabcur_north(2:n-1, 3:ndab,:)) ...
        - s*dt^2*udabcur_north(2:n-1, 2:ndab-1,:);

Next, we copy in the topmost row of points that the interior was able to update
into the first level of the auxiliary variables. ::
    
    % copy data to DAB boundary for
    % the top boundary in the y-direction
    udabnew_north(:,1,1) = unew(:,n-1);

Now, we run the CRBC recursions in the increasing direction of the auxiliary 
index to get updates to the bottommost points in the DAB layer. ::
    
    % run the "forward" recursion --- from eqn. 60-61 (a=ab=1,sig=sigb=0)
    w = 1/dt + c/h;
    
    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew_north(2:n-1,1,q+1) = ...
            (ab(q) - c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_north(2:n-1,1,q+1) ...
            +(ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_north(2:n-1,2,q+1) ...
            +(-a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_north(2:n-1,2,q) ...
            +(-a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_north(2:n-1,1,q) ...
            +(-ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_north(2:n-1,2,q+1) ...
            +(a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_north(2:n-1,2,q) ...
            +(a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_north(2:n-1,1,q);
    end

We begin the CRBC recursions at the topmost points at the highest auxilliary order by
applying the Sommerfeld radiation condition. Then we run the CRBC recursions in
decreasing auxiliary order. ::
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew_north(2:n-1,ndab, p+1) = ((udabcur_north(2:n-1,ndab-1, p+1) - udabnew_north(2:n-1,ndab-1, p+1)...
        + udabcur_north(2:n-1,ndab, p+1)) / dt + c*(udabcur_north(2:n-1,ndab-1,p+1) ...
        - udabcur_north(2:n-1,ndab, p+1) + udabnew_north(2:n-1,ndab-1, p+1))/h)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        udabnew_north(2:n-1,ndab,q) = ...
            (a(q) - c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_north(2:n-1,ndab,q) ...
            +(a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_north(2:n-1,ndab-1,q) ...
            +(-ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_north(2:n-1,ndab-1,q+1) ...
            +(-ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_north(2:n-1,ndab,q+1) ...
            +(-a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabnew_north(2:n-1,ndab-1,q) ...
            +(ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_north(2:n-1,ndab-1,q+1) ...
            +(ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_north(2:n-1,ndab,q+1);
    end

Next, we update the DAB layer on the right side of the domain. Note that these
updates are essentially identical to the updates on the top, but the indices 
are swapped because the recursions at the top have y-derivatives and the 
recursions on the right have x-derivatives (corresponding to the outward normal
directions at the respective boundary). It is possible to use a different number
of recursions on each face, but we do not do that here for simplicity.

To begin the DAB update, 
in all auxiliary variables we use the same update equation that we use to evolve
the interior points. ::
    
    %% East DAB
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    udabnew_east(2:ndab-1, 2:n-1,:) = 2*udabcur_east(2:ndab-1, 2:n-1,:) - udabold_east(2:ndab-1, 2:n-1, :) ...
        + ((c*dt)/h)^2 * (udabcur_east(3:ndab, 2:n-1, :) - 4*udabcur_east(2:ndab-1, 2:n-1,:) ...
        + udabcur_east(1:ndab-2, 2:n-1,:) + udabcur_east(2:ndab-1, 1:n-2,:) + udabcur_east(2:ndab-1, 3:n,:))...
        - s*dt^2*udabcur_east(2:ndab-1, 2:n-1,:);
    

Next, we copy in the rightmost column of points that the interior was able to update
into the first level of the auxiliary variables. ::
    
    % copy data to DAB boundary for
    % the right boundary in the x-direction
    udabnew_east(1,:,1) = unew(n-1,:);

Now, we run the CRBC recursions in the increasing direction of the auxiliary 
index to get updates to the leftmost points in the DAB layer. ::
    
    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew_east(1,2:n-1,q+1) = ...
            (ab(q) - c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_east(1,2:n-1,q+1) ...
            +(ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_east(2,2:n-1,q+1) ...
            +(-a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_east(2,2:n-1,q) ...
            +(-a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_east(1,2:n-1,q) ...
            +(-ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_east(2,2:n-1,q+1) ...
            +(a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_east(2,2:n-1,q) ...
            +(a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_east(1,2:n-1,q);
    end

We begin the CRBC recursions at the leftmost points at the highest auxilliary order by
applying the Sommerfeld radiation condition. Then we run the CRBC recursions in
decreasing auxiliary order. ::
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew_east(ndab,2:n-1, p+1) = ((udabcur_east(ndab-1,2:n-1, p+1) - udabnew_east(ndab-1,2:n-1, p+1)...
        + udabcur_east(ndab,2:n-1, p+1)) / dt + c*(udabcur_east(ndab-1,2:n-1,p+1) ...
        - udabcur_east(ndab,2:n-1, p+1) + udabnew_east(ndab-1,2:n-1, p+1))/h)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        udabnew_east(ndab,2:n-1,q) = ...
            (a(q) - c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_east(ndab,2:n-1,q) ...
            +(a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_east(ndab-1,2:n-1,q) ...
            +(-ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_east(ndab-1,2:n-1,q+1) ...
            +(-ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_east(ndab,2:n-1,q+1) ...
            +(-a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabnew_east(ndab-1,2:n-1,q) ...
            +(ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_east(ndab-1,2:n-1,q+1) ...
            +(ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_east(ndab,2:n-1,q+1);
    end

Next, we perform the updates for the corner DAB layer. The idea here is that we
have a doubly indexed set of auxiliary variables, where one auxiliary index 
corresponds to recursions running in the x-direction and the second index 
corresponds to recursions in the y-direction.

We begin the updates as before by applying the internal update equation to all of
the auxilliary vairables. ::

  %% Northeast Corner DAB
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    udabnew_ne(2:ndab-1, 2:ndab-1,:,:) = 2*udabcur_ne(2:ndab-1, 2:ndab-1,:,:) - udabold_ne(2:ndab-1, 2:ndab-1, :,:) ...
        + ((c*dt)/h)^2 * (udabcur_ne(3:ndab, 2:ndab-1, :,:) - 4*udabcur_ne(2:ndab-1, 2:ndab-1,:,:) ...
        + udabcur_ne(1:ndab-2, 2:ndab-1,:,:) + udabcur_ne(2:ndab-1, 1:ndab-2,:,:) + udabcur_ne(2:ndab-1, 3:ndab,:,:))...
        - s*dt^2*udabcur_ne(2:ndab-1, 2:ndab-1,:,:);

Next, we copy in the rightmost layer of auxiliary variables the top DAB layer was
able to update and the topmost layer of auxiliary variables the right DAB layer
was able to update. ::

    % copy data to DAB boundary from the North and East DAB layers
    udabnew_ne(2:ndab-1, 1,:,1) = udabnew_east(2:ndab-1, n-1, :); % first aux index corresponds to x-direction
    udabnew_ne(1, 2:ndab-1,1,:) = udabnew_north(n-1, 2:ndab-1, :); % second aux index corresponds to y-direction

Then, we run the forward recursions. Note that we have to run the recursions in 
the y-direction to get updates for the bottom row of points in the DAB layer and
recursions in the x-direction to get the updates for the points on the left side
of the layer. :: 

  % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        
        % recursions in the x-direction
        udabnew_ne(1,2:ndab-1,q+1,:) = ...
            (ab(q) - c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(1,2:ndab-1,q+1,:) ...
            +(ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(2,2:ndab-1,q+1,:) ...
            +(-a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(2,2:ndab-1,q,:) ...
            +(-a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(1,2:ndab-1,q,:) ...
            +(-ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_ne(2,2:ndab-1,q+1,:) ...
            +(a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_ne(2,2:ndab-1,q,:) ...
            +(a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_ne(1,2:ndab-1,q,:);
        
        % recursions in the y-direction
        udabnew_ne(2:ndab-1,1,:,q+1) = ...
            (ab(q) - c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(2:ndab-1,1,:,q+1) ...
            +(ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(2:ndab-1,2,:,q+1) ...
            +(-a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(2:ndab-1,2,:,q) ...
            +(-a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabcur_ne(2:ndab-1,1,:,q) ...
            +(-ab(q) + c*dt/h - sigb(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_ne(2:ndab-1,2,:,q+1) ...
            +(a(q) + c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_ne(2:ndab-1,2,:,q) ...
            +(a(q) - c*dt/h + sig(q))/(ab(q) + c*dt/h + sigb(q)) * udabnew_ne(2:ndab-1,1,:,q);
    end

Next, we apply the Sommerfeld termination condition and run the backward recurion 
updates. Again, note that to get the updates at the top of the layer, we use the
termination and recursions in the y-direction and to the updates at the right
we use the terminations and recursions in the x-direction. ::

  % apply the termination conditon, sommerfeld --- from eqn 56-57
    % x-direction
    udabnew_ne(ndab,2:ndab-1,p+1,:) = ((udabcur_ne(ndab-1,2:ndab-1,p+1,:) - udabnew_ne(ndab-1,2:ndab-1,p+1,:)...
        + udabcur_ne(ndab,2:ndab-1,p+1,:)) / dt + c*(udabcur_ne(ndab-1,2:ndab-1,p+1,:) ...
        - udabcur_ne(ndab,2:ndab-1,p+1,:) + udabnew_ne(ndab-1,2:ndab-1,p+1,:))/h)/w;
    
    % y-direction
    udabnew_ne(2:ndab-1,ndab,:,p+1) = ((udabcur_ne(2:ndab-1,ndab-1,:,p+1) - udabnew_ne(2:ndab-1,ndab-1,:,p+1)...
        + udabcur_ne(2:ndab-1,ndab,:,p+1)) / dt + c*(udabcur_ne(2:ndab-1,ndab,:,p+1) ...
        - udabcur_ne(2:ndab-1,ndab,:,p+1) + udabnew_ne(2:ndab-1,ndab-1,:,p+1))/h)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        
        % recursions in the x-direction
        udabnew_ne(ndab,2:ndab-1,q,:) = ...
            (a(q) - c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(ndab,2:ndab-1,q,:) ...
            +(a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(ndab-1,2:ndab-1,q,:) ...
            +(-ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(ndab-1,2:ndab-1,q+1,:) ...
            +(-ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(ndab,2:ndab-1,q+1,:) ...
            +(-a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabnew_ne(ndab-1,2:ndab-1,q,:) ...
            +(ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_ne(ndab-1,2:ndab-1,q+1,:) ...
            +(ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_ne(ndab,2:ndab-1,q+1,:);
        
        % recursions in the y-direction
        udabnew_ne(2:ndab-1,ndab,:,q) = ...
            (a(q) - c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(2:ndab-1,ndab,:,q) ...
            +(a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(2:ndab-1,ndab-1,:,q) ...
            +(-ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(2:ndab-1,ndab-1,:,q+1) ...
            +(-ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabcur_ne(2:ndab-1,ndab,:,q+1) ...
            +(-a(q) + c*dt/h - sig(q))/(a(q) + c*dt/h + sig(q)) * udabnew_ne(2:ndab-1,ndab-1,:,q) ...
            +(ab(q) + c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_ne(2:ndab-1,ndab-1,:,q+1) ...
            +(ab(q) - c*dt/h + sigb(q))/(a(q) + c*dt/h + sig(q)) * udabnew_ne(2:ndab-1,ndab,:,q+1);
    end


Finally, we copy the updated auxiliary variables from the corner to the two sides
and then we copy the updated first level auxiliary variables from the sides
into the internal solver. ::
    
    % copy values from the Northest DAB corner to the North and East DAB
    % layers
    udabnew_east(2:ndab-1, n, :) = udabnew_ne(2:ndab-1, 2,:,1);
    udabnew_north(n, 2:ndab-1, :) = udabnew_ne(2, 2:ndab-1,1,:);
    
    % copy values from the DAB layers into the internal solver
    unew(n,:) = udabnew_east(2,:,1);
    unew(:,n) = udabnew_north(:,2,1);

Notice that we never update the corners of the corner DAB layer. We do this because
those values are never used in any of the updates, but it is also fortuitous
because the corner points are overdetermined and it is not always clear how to 
deal with this.
    
We plot the field values and the error by comparing to the larger simulation. :: 

    % figures
    
    % field and comparison to larger simulation
    if (mod(t,plot_freq) == 0)
        subplot(1,2,1)
        surf(unew'); % transpose because matlab plots x and y "backwards"
        view(2)
        shading interp
        colorbar;
        title('field values')
        subplot(1,2,2)
        surf(unew' - runew(1:n,1:n)')
        colorbar;
        view(2)
        shading interp
        title('Error compared to larger simulation')
        drawnow
    end


Lastly, we rotate the storage arrays so we can procede to the next time step. ::
    
    % swap old, new, and current values
    uold = ucur;
    ucur = unew;
    
    ruold = rucur;
    rucur = runew;
    
    udabold_east = udabcur_east;
    udabcur_east = udabnew_east;
    
    udabold_north = udabcur_north;
    udabcur_north = udabnew_north;
    
    udabold_ne = udabcur_ne;
    udabcur_ne = udabnew_ne;
    
  end


.. rubric:: References

.. bibliography:: zcite.bib
   :encoding: UTF8
   :list: enumerated
   :filter: author % "Givoli"
