.. highlight:: matlab

***************
1D Klein-Gordon
***************

This tutorial explains the code `klein_gordon_1d_DAB.m <https://bitbucket.org/rbcpack/rbcpack/src/default/YeeCRBC/Matlab/klein_gordon_1d_DAB.m>`_ and requires the use of `optimal_cosines.m <https://bitbucket.org/rbcpack/rbcpack/src/default/YeeCRBC/Matlab/optimal_cosines.m>`_

Introduction
============

This Matlab code implements a second order finite difference approximation to
the 1D Klien-Gordon equation. On one side, the grid is terminated with a Double
Absorbing Boundary (DAB).

We use the Klein-Gordon equation instead of the wave equation because the 
Sommerfeld radiation condition is the correct boundary condition for the wave
equation in 1D. For Klein-Gordon, dispersion results in waves moving at different
speeds, which means Sommerfeld is no longer exact and the benefits of the DAB can
be better illustrated.


What this program does
----------------------

For this example, we consider the 1D Klein-Gordon equation,

.. math::
  :nowrap:

  \begin{align}
    \frac{d^2 u}{d t^2} &= c^2 \frac{d^2 u}{d x^2} - s u,
  \end{align}

where :math:`c>0` and :math:`s>0`.

For this example, we will impose Dirichlet boundary conditions on the left and
truncate the domain on the right with a DAB. We will use zero initial conditions.

.. math::
  :nowrap:

  \begin{align}
    u(x,0) = 0, \\
    \frac{du}{dt}(x,0) = 0, \\
    u(0,t) = 0,
  \end{align}

.. _discretization:

We discretize this equation using second order finite differences on the 
domain :math:`[0, 1]` with mesh spacing of :math:`\Delta x` and we define

.. math::
  :nowrap:

  \begin{align}
    x_i &=  i\Delta x, \qquad i=0,...,n-1.
  \end{align}

We choose a time step size, :math:`\Delta t`, satisfying

.. math::
  :nowrap:

  \begin{align}
    \Delta t \leq \frac{dx}{c}.
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
    \frac{u(x_i,t_{n-1}) - 2u(x_i,t_{n}) + u(x_i,t_{n+1})}{\Delta t} 
    = c\left( \frac{u(x_{i-1},t_{n}) - 2u(x_i,t_{n}) + u(x_{i+1},t_{n})}{\Delta x} \right)
    - s u(x_i,t_{n})
  \end{align}

We use the discretization of the DAB described in :ref:`theory overview <dab_discretization>`


We drive the simulation with a point source that takes the form of
a differentiated Gaussian.


The commented program
=====================

We begin by choosing some basic simulation parameters. First we choose the number
of grid points to use in the discetization. Then we choose the number of grid
points to extend the domain by so we can compare to a larger simulation to check
error. In this case, we use 500 grid points and extend the larger simulation by
600 grid points, which corresponds to running on the domain [0,2.2]. Then, we
choose the problem and source parameters. ::

  % domain parameters
  n = 500;           % number of grid points in domain [0,1]
  m = 600;           % number of grid points to extent the domain by for a reference 
                     % solution using a larger simulation
  s = 25;            % dispersion parameter
  c = 1;             % wave speed
  nsteps = 1500;     % number of time steps
  cfl = 0.99;        % cfl ratio 

  % compute grid spacing / time step
  dx = 1.0 / (n - 1);
  dt = cfl * dx / c;

  % source paramters
  tw = 25*dt;             % pulse width
  t0 = 5*tw;              % pulse delay (roughly need t0 > 4*tw to be smooth)
  amp = 1;                % pulse "amplitude"
  sloc = 90;              % source location

Next, we choose the DAB parameters. We choose the number of recursions to use,
:math:`p` and how wide the DAB layer should be in grid points. We require at 
least three points to support the update to the Klein-Gordon equation. For
efficiency and accuracy, 3 is the best choice; however, if the auxiliary variables
are to be plotted, increasing the thickness is desireable.

We present two ways to choose the parameter values, the first is simplistic and
sets all of the values to either 1 or 0. This generally works well for short 
times. The alternative is to choose the optimal parameters, which provide an
error estimate valid until the provided time. ::

  % DAB parameters
  p = 3;                  % DAB/CRBC order
  ndab = 3;               % DAB width
  % a = ones(p,1);          % choose all the cosines to be one for simplicity
  % ab = ones(p,1);         % choose all the cosines to be one for simplicity
  % sig = zeros(p,1);       % choose all the sigmas to be zero for simplicity
  % sigb = zeros(p,1);      % choose all the sigmas to be zero for simplicity

  % ... or use optimal cosines
  T = nsteps * dt;
  eta = (n - sloc)*dx / (c * T);
  if (p>0)
    [at errest] = optimal_cosines(eta, p-1);
    a = at(1:2:2*p);
    ab = at(2:2:2*p);
    sig = 0.5*dt*(1 - a.*a) ./ (T*a);
    sigb = 0.5*dt*(1 - ab.*ab) ./ (T*ab);
  end

Next, we allocate the storage for all of the field values and auxiliary variables
we will use. ::

  % allocate storage
  unew = zeros(n,1); % field values
  ucur = zeros(n,1);
  uold = zeros(n,1);

  runew = zeros(n+m,1); % for larger reference simulation 
  rucur = zeros(n+m,1);
  ruold = zeros(n+m,1);

  udabnew = zeros(ndab, p+1); % dab aux. variables
  udabcur = zeros(ndab, p+1);
  udabold = zeros(ndab, p+1);

We begin time stepping by updating all of the internal field values and adding a
source term. ::

  % time step
  for t=1:nsteps
    
    % internal updates --- eqn 54, in DAB paper
    unew(2:n-1) = 2*ucur(2:n-1) - uold(2:n-1) + ((c*dt)/dx)^2 * (ucur(3:n) ...
        - 2*ucur(2:n-1) + ucur(1:n-2)) - s*dt^2*ucur(2:n-1);
    
    % reference solution
    runew(2:m+n-1) = 2*rucur(2:m+n-1) - ruold(2:m+n-1) + ((c*dt)/dx)^2 * (rucur(3:m+n) ...
        - 2*rucur(2:m+n-1) + rucur(1:m+n-2)) - s*dt^2*rucur(2:m+n-1);
    
    % add source 
    unew(sloc) = unew(sloc) - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);
    runew(sloc) = runew(sloc) - 2*((t*dt - t0)/tw)*amp*exp(-((t*dt - t0)/tw)^2);

To begin the DAB update, in all auxiliary variables we use the same update 
equation that we use to evolve the interior points. ::
    
    % perform wave equation update for the interior of the DAB --- eqn 54, in DAB paper
    for q=1:p+1
       udabnew(2:ndab-1, q) = 2*udabcur(2:ndab-1, q) - udabold(2:ndab-1, q) + ...
           ((c*dt)/dx)^2 * (udabcur(3:ndab, q)- 2*udabcur(2:ndab-1, q) +...
           udabcur(1:ndab-2, q)) - s*dt^2*udabcur(2:ndab-1, q); 
    end

Next, we copy in the rightmost point that the interior was able to update
into the first level of the auxiliary variables. ::
    
    % left boundary is correctly set to zero, copy data to DAB boundary for
    % the right boundary
    udabnew(1,1) = unew(n-1);

Now, we run the CRBC recursions in the increasing direction of the auxiliary 
index to get updates to the leftmost point in the DAB layer. ::
    
    % run the "forward" recursion --- from eqn. 60-61 (a=ab=1,sig=sigb=0)
    w = 1/dt + c/dx;

    % run the "forward" recursion --- from eqn. 60-61, generalized
    for q=1:p
        udabnew(1,q+1) = ...
            (ab(q) - c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(1,q+1) ...
            +(ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2,q+1) ...
            +(-a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(2,q) ...
            +(-a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabcur(1,q) ...
            +(-ab(q) + c*dt/dx - sigb(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2,q+1) ...
            +(a(q) + c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(2,q) ...
            +(a(q) - c*dt/dx + sig(q))/(ab(q) + c*dt/dx + sigb(q)) * udabnew(1,q);   
    end

We begin the CRBC recursions at the rightmost at the highest auxilliary order by
applying the Sommerfeld radiation condition. Then we run the CRBC recursions in
decreasing auxiliary order. ::
    
    % apply the termination conditon, sommerfeld --- from eqn 56-57
    udabnew(ndab, p+1) = ((udabcur(ndab-1, p+1) - udabnew(ndab-1, p+1) + udabcur(ndab, p+1)) / dt ...
        + c*(udabcur(ndab-1,p+1) - udabcur(ndab, p+1) + udabnew(ndab-1, p+1))/dx)/w;

    % run the "backward" recursions --- from eqn. 58-59, generalized\
    for q=p:-1:1
        udabnew(ndab,q) = ...
            (a(q) - c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab,q) ...
            +(a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab-1,q) ...
            +(-ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab-1,q+1) ...
            +(-ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabcur(ndab,q+1) ...
            +(-a(q) + c*dt/dx - sig(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(ndab-1,q) ...
            +(ab(q) + c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(ndab-1,q+1) ...
            +(ab(q) - c*dt/dx + sigb(q))/(a(q) + c*dt/dx + sig(q)) * udabnew(ndab,q+1);   
    end

Finally, we copy the updated first level auxiliary variable into the internal
solver. ::
    
    % copy DAB value back into the field
    unew(n) = udabnew(2,1);
    
We plot the field values and the error by comparing to the larger simulation. 
The commented out portion plots the field values and the auxiliary layers (these
plots are clearer if the DAB layer is relatively wide). ::

    % figures
    
    % field and comparison to larger simulation
    figure(1)
    subplot(1,2,1)
    plot(1:n, unew);
    title('field values')
    subplot(1,2,2)
    plot(1:n, unew - runew(1:n))
    title('Error compared to larger simulation')
    drawnow
    
    % field and auxiliary fields
  % figure(2)
  % subplot(1, p+4, 1:3)
  % plot(1:n, unew);
  % title('field values')
  % for i=1:p+1
  %   subplot(1, p+4, i+3)
  %   plot(1:ndab, udabnew(:,i));
  %   title(sprintf('p = %i', i-1))
  % end
  % drawnow

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
