.. highlight:: matlab

****************
2D Wave Equation 
****************

This tutorial explains the code `wave_2d_DAB.m <https://github.com/jrlagrone/rbcpack/blob/main/YeeCRBC/Matlab/optimal_cosines.m>`_

Introduction
============

This Matlab code implements a second order finite difference approximation to
the 2D wave equation. On one side, the grid is terminated with a Double
Absorbing Boundary (DAB).

The primary thing to notice here is that the DAB is essentially identical to the
1D case described in the :doc:`1D Klein-Gordon example <1d_kg_matlab>`. This is
because we only need to use the appropriate wave equation to update the interior
of the layer and all remaining updates are in the normal direction.


What this program does
----------------------

For this example, we consider the 2D wave equation,

.. math::
  :nowrap:

  \begin{align}
    \frac{d^2 u}{d t^2} &= c^2 \left(\frac{d^2 u}{d x^2} + \frac{d^2 u}{d y^2} \right),
  \end{align}

where :math:`c>0`.

For this example, we will impose Dirichlet boundary conditions on the both sides
in the x-direction and at the bottom in the y-direction. We truncate the domain 
at the top in the y-direction with a DAB. We will use zero initial conditions.

.. math::
  :nowrap:

  \begin{align}
    u(x,y,0) = 0, \\
    \frac{du}{dt}(x,y,0) = 0, \\
    u(0,y,t) = 0, \\
    u(0.1, y,t) = 0, \\
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
  \end{align}

We use the discretization of the DAB described in :ref:`theory overview <dab_discretization>`


We drive the simulation with a point source that takes the form of
a differentiated Gaussian.


The commented program
=====================

We begin by choosing some basic simulation parameters. First we choose the number
of grid points to use in the discetization. Then we choose the number of grid
points to extend the domain by so we can compare to a larger simulation to check
error. In this case, we use 300 grid points and extend the larger simulation by
300 grid points, which corresponds to running on the domain [0,0.2]. Then, we
choose the problem and source parameters. ::

  % domain parameters
  n = 300;           % number of grid points in domain [0,0.1]
  m = 300;           % number of grid points to extent the domain by for a reference 
                   % solution using a larger simulation
  c = 1;             % wave speed
  nsteps = 700;     % number of time steps
  cfl = 0.99; % cfl ratio (1 is exact dispersion relation, but num. unstable)

  % compute grid spacing / time step
  dx = 0.1 / (n - 1); % we'll just use dx=dy
  dt = cfl / (c*sqrt(2/dx^2));

  % source paramters
  tw = 25*dt;             % pulse width
  t0 = 5*tw;              % pulse delay (roughly need t0 > 4*tw to be smooth)
  amp = 1;                % pulse "amplitude"
  sloc = [210, 175];      % source location

Next, we choose the DAB parameters. We choose the number of recursions to use,
:math:`p` and how wide the DAB layer should be in grid points. We require at 
least three points to support the update to the wave equation. For
efficiency and accuracy, 3 is the best choice; however, if the auxiliary variables
are to be plotted, increasing the thickness is desireable.

We present two ways to choose the parameter values, the first is simplistic and
sets all of the values to either 1 or 0. This generally works well for short 
times. The alternative is to choose the optimal parameters, which provide an
error estimate valid until the provided time. ::

  % DAB parameters
  p = 5;                  % DAB/CRBC order
  ndab = 3;               % DAB width
  a = ones(p,1);          % choose all the cosines to be one for simplicity
  ab = ones(p,1);         % choose all the cosines to be one for simplicity
  sig = zeros(p,1);       % choose all the sigmas to be zero for simplicity
  sigb = zeros(p,1);      % choose all the sigmas to be zero for simplicity

  % or choose optimal cosines
  T = nsteps * dt;
  eta = (n - sloc(2))*dx / (c * T);
  if (p>0)
    [at, errest] = optimal_cosines(eta, p-1);
    a = at(1:2:2*p);
    ab = at(2:2:2*p);
    sig = 0.5*dt*(1 - a.*a) ./ (T*a);
    sigb = 0.5*dt*(1 - ab.*ab) ./ (T*ab);
  end

Next, we allocate the storage for all of the field values and auxiliary variables
we will use. ::

  % allocate storage
  unew = zeros(n); % field values
  ucur = zeros(n);
  uold = zeros(n);

  runew = zeros(n, n+m); % for larger reference simulation 
  rucur = zeros(n, n+m);
  ruold = zeros(n, n+m);

  udabnew = zeros(n, ndab, p+1); % dab aux. variables
  udabcur = zeros(n, ndab, p+1);
  udabold = zeros(n, ndab, p+1);

We begin time stepping by updating all of the internal field values and adding a
source term. ::

  % time step
  for t=1:nsteps
    
    % internal updates --- eqn 54, in DAB paper
    unew(2:n-1, 2:n-1) = 2*ucur(2:n-1, 2:n-1) - uold(2:n-1, 2:n-1) ...
        + ((c*dt)/dx)^2 * (ucur(3:n, 2:n-1) - 4*ucur(2:n-1,2:n-1) ...
        + ucur(1:n-2, 2:n-1) + ucur(2:n-1, 1:n-2) + ucur(2:n-1, 3:n));
    
    % reference solution (same thing on a bigger domain)
    runew(2:n-1, 2:n+m-1) = 2*rucur(2:n-1, 2:n+m-1) - ruold(2:n-1, 2:n+m-1) ...
        + ((c*dt)/dx)^2 * (rucur(3:n, 2:n+m-1) - 4*rucur(2:n-1,2:n+m-1) ...
        + rucur(1:n-2, 2:n+m-1) + rucur(2:n-1, 1:n+m-2) + rucur(2:n-1, 3:n+m));
    
    % add source
    unew(sloc(1), sloc(2)) = unew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    runew(sloc(1), sloc(2)) = runew(sloc(1), sloc(2)) ...
        - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);

To begin the DAB update, in all auxiliary variables we use the same update 
equation that we use to evolve the interior points. ::
    
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    udabnew(2:n-1, 2:ndab-1,:) = 2*udabcur(2:n-1, 2:ndab-1,:) - udabold(2:n-1, 2:ndab-1,:) ...
        + ((c*dt)/dx)^2 * (udabcur(3:n, 2:ndab-1,:) - 4*udabcur(2:n-1,2:ndab-1,:) ...
        + udabcur(1:n-2, 2:ndab-1,:) + udabcur(2:n-1, 1:ndab-2,:) + udabcur(2:n-1, 3:ndab,:));

Next, we copy in the topmost row of points that the interior was able to update
into the first level of the auxiliary variables. ::
    
    % copy data to DAB boundary for
    % the right boundary in the y-direction
    udabnew(:,1,1) = unew(:,n-1);

Now, we run the CRBC recursions in the increasing direction of the auxiliary 
index to get updates to the bottommost points in the DAB layer. ::
    
    % run the "forward" recursion --- from eqn. 60-61 (a=ab=1,sig=sigb=0)
    w = 1/dt + c/dx;

    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew(2:n-1,1,q+1) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,1,q+1) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,2,q+1) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,2,q) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2:n-1,1,q) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2:n-1,2,q+1) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2:n-1,2,q) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2:n-1,1,q);   
    end

We begin the CRBC recursions at the topmost points at the highest auxilliary order by
applying the Sommerfeld radiation condition. Then we run the CRBC recursions in
decreasing auxiliary order. ::
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew(2:n-1,ndab, p+1) = ((udabcur(2:n-1,ndab-1, p+1) - udabnew(2:n-1,ndab-1, p+1) + udabcur(2:n-1,ndab, p+1)) / dt ...
        + c*(udabcur(2:n-1,ndab-1,p+1) - udabcur(2:n-1,ndab, p+1) + udabnew(2:n-1,ndab-1, p+1))/dx)/w;
    
    % run the "backward" recursions --- from eqn. 58-59, generalized
    for q=p:-1:1
        udabnew(2:n-1,ndab,q) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab,q) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab-1,q) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab-1,q+1) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(2:n-1,ndab,q+1) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(2:n-1,ndab-1,q) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(2:n-1,ndab-1,q+1) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(2:n-1,ndab,q+1);   
    end

Finally, we copy the updated first level auxiliary variables into the internal
solver. ::
    
    % copy DAB value back into the field
    unew(:,n) = udabnew(:,2,1);
    
We plot the field values and the error by comparing to the larger simulation. 
The commented out portion plots the field values and the auxiliary layers (these
plots are clearer if the DAB layer is relatively wide). ::

    % figures
    
    % field and comparison to larger simulation
    figure(1)
    subplot(1,2,1)
    surf(unew);
    view(2)
    shading interp
    colorbar;
    title('field values')
    subplot(1,2,2)
    surf(unew - runew(1:n,1:n))
    colorbar;
    view(2)
    shading interp
    title('Error compared to larger simulation')
    drawnow
    
    % field and auxiliary fields
  %     figure(2)
  %     subplot(1, p+4, 1:3)
  %     surf(unew);
  %     view(2)
  %     shading interp
  %     colorbar;
  %     title('field values')
  %     for i=1:p+1
  %         subplot(1, p+4, i+3)
  %         surf(udabnew(:,:,i));
  %         shading interp
  %         view(2)
  %         colorbar;
  %         title(sprintf('p = %i', i-1))
  %     end
  %     drawnow   

Lastly, we rotate the storage arrays so we can procede to the next time step. ::
    
    % swap old, new, and current values
    uold = ucur;
    ucur = unew;
    
    ruold = rucur;
    rucur = runew;
    
    udabold = udabcur;
    udabcur = udabnew;
    
  end


.. rubric:: References

.. bibliography:: zcite.bib
   :encoding: UTF8
   :list: enumerated
   :filter: author % "Givoli"
