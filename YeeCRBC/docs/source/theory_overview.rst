**************************
A Brief Overview of Theory
**************************

Overview of Complete Radiation Boundary Conditions
==================================================

Following from :cite:`CRBC`, the goal is first write down a complete
representation of the solution away from the interior of our computational
domain.
To this end, suppose that  :math:`u(x,y,t)` satisfies 

.. math::
  :label: half_space_we
  
      &\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u, \quad x > - \delta, 
        \quad y \in \mathbb{R}^{d-1}, \quad t > 0, \\
      &u(x,y,0) = 0,
  
for some  :math:`\delta > 0`.
Furthermore, suppose that the field is produced sources, scatterers, and other inhomogeneities located
in the half space  :math:`x < -\delta` and we can account for these effects with the Dirichlet data:

.. math::
  :label: dir_data

    u(-\delta, y, t) = g(y,t). 

Taking the transverse-Fourier-Laplace transforms, we can write the solution to 
the problem :eq:`half_space_we`-:eq:`dir_data` as

.. math::

  \hat{u}(x,k,s) = \hat{u}(-\delta,k,s) e^{-(\bar{s}^2 + |k|^2)^{1/2}(x+\delta)},

where  :math:`k` are the dual Fourier variables to the transverse spatial coordinates,
:math:`| \cdot |` is the Euclidean norm,  :math:`s` is the dual Laplace variable
to time, and :math:`\bar{s} = \frac{s}{c}`.
The branch of the square root is chosen to have positive real part when the real 
part of  :math:`s` is positive.
We invert the transformed solution with respect to  :math:`k` over the hyperplane  
:math:`\mathbb{R}^{d-1}` and with respect to  :math:`s` over the contour  
:math:`s = i \omega + \frac{1}{T}` with  :math:`\omega` real and fixed  :math:`T`
(with units of time).
To accomplish this, set

.. math::

  k= \frac{\tilde{k}}{cT}, \qquad \omega = \frac{\tilde{\omega}}{T},

then

.. math::
  :label: sqrt

    (\bar{s}^2 + |k|^2)^{1/2} = ((1+i\tilde{\omega})^2 + \tilde{k}^2)^{1/2} = a + ib,

for some  :math:`a,b \in \mathbb{R}`. Squaring both sides, we get

.. math::

  b = \frac{\tilde{\omega}}{a}, \qquad 1 + \tilde{k}^2 - \tilde{\omega}^2 = a^2 - \frac{\tilde{\omega^2}}{a^2}.

In order to agree with the chosen branch of the square root, we require
:math:`a \geq 1`. Therefore, for some  :math:`\phi \in [0, \frac{\pi}{2})`, 
we can write

.. math::

  a= \frac{1}{\cos \phi}, \qquad b = \tilde{\omega} \cos \phi.

Substituting back into :eq:`sqrt`,

.. math::
  :nowrap:
  
  \begin{align}
    (\bar{s}^2 + |k|^2)^{1/2} = \bar{s} \cos \phi + \frac{1}{cT} \frac{\sin^2 \phi}{\cos \phi}.
  \end{align} 

Finally, writing the inverse transforms in polar coordinates, 
:math:`\rho, \theta \in \mathbb{S}^{d-2}`, in the dual space variables and 
replacing  :math:`\rho = |k|` by  :math:`\rho(\phi, \omega)`, we have

.. math::
  :nowrap:
  
  \begin{align}
    u(x,y,t) = (2 \pi)^{-\frac{d+1}{2}} \int_{0}^{\frac{\pi}{2}} \int_{-\infty}^{\infty}
               \int_{\mathbb{S}^{d-2}} e^{\psi} \bar{u} \rho^{d-2} \frac{\partial \rho}{\partial \phi}
                dA(\theta) d\omega d\phi,
  \end{align} 

where

.. math::
  :nowrap: 
  
  \begin{align}
     \psi &= \left(i\omega + \frac{1}{T} \right) \left(t- \frac{\cos \phi}{c}(x+\delta) \right)
            + i \rho(\theta \cdot y) - \frac{1}{cT}\frac{\sin^2 \phi}{\cos \phi}(x+\delta), \\
     \bar{u} &= \hat{u}\left(-\delta, \rho \theta, \frac{1}{T+i \omega} \right).
  \end{align} 

Setting

.. math::
  :nowrap:
  
  \begin{align}
    \Phi(t,y,\phi) =  (2 \pi)^{-\frac{d+1}{2}}  \int_{-\infty}^{\infty} \int_{\mathbb{S}^{d-2}} 
                      e^{\left(i\omega + \frac{1}{T} \right) + i \rho(\theta \cdot y)} \bar{u} \rho^{d-2}
                      \frac{\partial \rho}{\partial \phi} dA(\theta) d\omega,
  \end{align} 

we can write the complete wave representation of the solution, valid for  :math:`x > -\delta`,
as

.. math::
  :label: complete_wave_rep
  
    u(x,y,t) = \int_{0}^{\frac{\pi}{2}} \Phi \left( t - \frac{\cos \phi}{c}(x+\delta), y, \phi \right)
               e^{- \frac{1}{cT}\frac{\sin^2 \phi}{\cos \phi}(x+\delta)} d \phi.

To arrive at an approximate local boundary condition, we can approximate the  :math:`\phi`
integral in :eq:`complete_wave_rep` by an appropriate quadrature rule using
nodes  :math:`\phi_j` and weights  :math:`h_j`:
  
.. math::
  :label: approx_complete_wave_rep
  
  u(x,y,t) \approx \sum\limits_{j=0}^{P} h_j \Phi \left( t - \frac{\cos \phi_j}{c}(x+\delta), y, \phi_j \right) 
            e^{- \frac{1}{cT}\frac{\sin^2 \phi_j}{\cos \phi_j}(x+\delta)}.


Introducing a second set of angles  :math:`\bar{\phi}_j` and defining auxiliary functions 
:math:`u_j(x,y,t)` by setting  :math:`u_0 = u` and recursively solving

.. math::
  :label: crbc_recursion

  \bar{a}_{j} \frac{\partial u_{j+1}}{\partial t} - \frac{\partial u_{j+1}}{\partial x} + \bar{\sigma}_{j} 
     u_{j+1}  = a_j \frac{\partial u_{j}}{\partial t} + \frac{\partial u_{j}}{\partial x} + \sigma_{j} u_{j}, 

where

.. math::
  :nowrap:
  
  \begin{align}
    a_j = \frac{cos \phi_j}{c}, \qquad \bar{a}_j = \frac{cos \bar{\phi}_j}{c}, \qquad \sigma_j 
        = \frac{1}{cT}\frac{\sin^2 \phi_j}{\cos \phi_j}, \qquad \bar{\sigma}_j 
        = \frac{1}{cT}\frac{\sin^2 \bar{\phi}_j}{\cos \bar{\phi}_j}, 
  \end{align} 

subject to  :math:`u_{j+1}(x,y,0) = 0`. Substituting :eq:`approx_complete_wave_rep` for  :math:`u_0`, 
we notice that individual terms of :eq:`approx_complete_wave_rep` are annihilated by the right-hand 
side of :eq:`crbc_recursion`; therefore, 

.. math::
  :label: crbc_p

    u_{P+1} = 0.

Finally, to utilize this approximation, we impose :eq:`crbc_p` on the incoming
normal characteristic variables.
  
Advantages
----------

The complete radiation boundary conditions (CRBCs) have several advantages over
other boundary conditions

* The boundaries are local, but achieve comparable long-time accuracies to nonlocal conditions. This allows for cheaper computations and a greater range of geometries to be considered.

* The number of auxiliary variables,  :math:`P`, required to obtain an accuracy  :math:`\varepsilon` up to time  :math:`T` satisfies  :math:`P = O\left(\ln \frac{1}{\varepsilon} \cdot \ln \frac{cT}{\delta}\right)`. This gives a clear notion of accuracy as well as providing a means to select the optimal parameters `a priori`.

* The recursions :eq:`crbc_recursion` require only first derivatives even though :math:`P` may be arbitrarily high, which enables the boundary conditions to be easily implemented to an arbitrary order :math:`P` for a given discretization scheme.

Overview of the Double Absorbing Boundary
=========================================

Following from :cite:`DAB`, we illustrate the double absorbing boundary layer (DAB)
using the Klein-Gordon equation in the semi-infinite wave guide:

.. math::
  :nowrap: 
  
  \begin{align}
     Wu \equiv \frac{\partial^2 u}{\partial t^2} - c^2 \nabla^2 u + su = f, \\
     u(x_L,y,t) = u(x,y_L,t) = u(x,y_R,t) = 0, \\
     u(x,y,0) = g(x,y), \\
     \frac{\partial u}{\partial t} (x,y,0) = \dot{g}(x,y). 
  \end{align} 

We further suppose that if :math:` x > x_I > x_L`, the medium is homogeneous and free of sources, that is,
`c` and :math:`s` are constants and :math:`f = 0`; therefore, :math:`Wu = 0`. 
Additionally, we require that the initial conditions vanish so :math:`g = \dot{g} = 0`
for :math:`x > x_I`. 

We now truncate the domain at some :math:`x = x_R > x_I`. We can view the entire truncated domain
as being divided into two sub-domains: the interior domain :math:`\varOmega_I \equiv [x_L, x_I] \times [y_L, y_R]` and a thin layer :math:`\varOmega_L \equiv [x_I, x_R] \times [y_L, y_R]`. The goal is to 
construct a solution in :math:`\varOmega_I` as close as possible to the solution of the original semi-infinite 
problem in that domain and to use :math:`\varOmega_L` as an absorbing layer.

To do this, we introduce auxiliary variables :math:`u_0, ... u_{P+1}` in :math:`\varOmega_L` and
require :math:`u_j` to satisfy that same wave equation as :math:`u`:

.. math::
  :nowrap: 
  
  \begin{align}
    Wu_j \equiv \frac{\partial^2 u_j}{\partial t^2} - c^2 \nabla^2 u_j + su_j = 0, \qquad j=0,1...,P+1, 
         \qquad \text{in }   \varOmega_L.
  \end{align} 

The auxiliary variables satisfy zero initial conditions

.. math::
  :nowrap: 
  
  \begin{align}
    u_j(x,y,0) = \frac{\partial u_j}{\partial t} =0, \qquad j=0,1...,P+1, \qquad \text{in } \varOmega_L, 
  \end{align} 

and satisfy the boundary conditions

.. math::
  :nowrap: 
  
  \begin{align}
    u_j(x,y,t) = 0, \qquad j=0,1...,P+1, \qquad y=y_L,y_R, \qquad \text{in } \varOmega_L, 
  \end{align} 

To define the additional boundary conditions on :math:`\varOmega_L`, we utilize the CRBC boundary recursions :eq:`crbc_recursion`:

.. math::
  :nowrap: 
  
  \begin{align}
     \bar{a}_{j} \frac{\partial u_{j+1}}{\partial t} - \frac{\partial u_{j+1}}{\partial x} + \bar{\sigma}_{j} u_{j+1}
      = a_j \frac{\partial u_{j}}{\partial t} + \frac{\partial u_{j}}{\partial x} + \sigma_{j} u_{j},
  \end{align} 

and require them to hold at :math:`x=x_I,x_R` (note that in principle, we could use any radiation/absorbing boundary conditions here). 
On :math:`x = x_I`, we require the :math:`u` and :math:`u_0` to agree in value
and slope:

.. math::
  :nowrap: 
  
  \begin{align}
    u_0 = u, \qquad \frac{\partial u_{0}}{\partial x} = \frac{\partial u}{\partial x}, \qquad x = x_I.
  \end{align} 

Since :math:`u` and :math:`u_0` satisfy the same wave equation in :math:`\varOmega_L`,

.. math::
  :nowrap: 
  
  \begin{align}
    u_0 \equiv u, \qquad \text{in } \varOmega_L.
  \end{align} 

Finally, at :math:`x = x_R`, we require the "termination condition"

.. math::
  :nowrap: 
  
  \begin{align}
    \bar{a}_{P+1} \frac{\partial u_{P+1}}{\partial t} - \frac{\partial u_{P+1}}{\partial x} 
    + \bar{\sigma}_{P+1} u_{P+1} = 0.
  \end{align} 


Advantages
----------

The following are some advantages the double absorbing boundary layer construction using CRBCs (note that is in principle to construct a DAB using any absorbing boundary conditions):

* The auxiliary variables are defined to satisfied the desired wave equation in the boundary layer.

* It is not necessary to eliminate normal derivatives in the formulation. This potentially makes the DAB easier to implement as apposed to directly using the CRBCs because we do not have to replace the normal derivatives with temporal and tangential derivatives. Furthermore, this may make it easier to implement in a layered medium.

* In some cases, it is not necessary to handle edges or corners in a special way.

* The interior of the DAB layer can be approximated using a different scheme than what is used on the domain :math:`\varOmega_I`, which allows us, for instance, to return data at whatever ghost nodes the scheme in :math:`\varOmega_I` may need.

Discretization of the DAB
=========================

The basic idea is to use the scalar wave equation to update the points in the
center a DAB layer and use CRBC recursions plus either a termination
condition or value from that can be updated by the solver in the interior of the
domain to update the outer points of the layer.
  
More precisely, supposing we have an initial condition compactly supported 
away from the boundary, to provide a boundary update for the left side in the
:math:`x` direction, we introduce the auxiliary variables 
:math:`u_{p}^{\tilde{i},j,k,n}` and :math:`u_{p}^{\tilde{i},j,k,n}` 
for :math:`p=0,...,P`  and :math:`\tilde{i} = n_x-1,n_x,n_x+1`.
The equations for the :math:`u` auxiliary variables are

.. math::

  u_{0} &= u, \qquad x =  x_{n_x-1} \\
  \frac{\partial^2 u_{p}}{\partial t^2} &= c^2 \nabla^2 u_{p}, \qquad x = x = x_{n_x}\\
  \left(a_p \frac{\partial}{\partial t} - c \frac{\partial }{\partial x} + \sigma_p  \right)u_{p-1}
  & = \left(\bar{a}_p \frac{\partial}{\partial t} + c \frac{\partial}{\partial x} + \bar{\sigma}_p \right) u_{p},
  \qquad x = x_{n_x-\frac{1}{2}}, x_{n_x+\frac{1}{2}},

with the termination condition

.. math::

  \left(\frac{\partial}{\partial t} -c \frac{\partial}{\partial x} \right)u_{P} = 0, & \quad x =  x_{n_x+\frac{1}{2}}.

We discretize these equations using second order centered differences. Additionally,
we average in space and time to ensure that the differences are centered at the
appropriate space-time coordinate. We note that this procedure is also second order.
This yields the following discretization:

.. math::

  u_{0}^{n_x-1,j,k,n} = E_{y}^{n_x-1,j,k,n}.

The discretization of the wave equation is

.. math::

  u_{p}^{n_x,j,k,n-1} - 2 u_{p}^{n_x,j,k,n} + u_{p}^{n_x,j,k,n+1}   & = 
  c^2 (\Delta t)^2 \left[ \frac{1}{(\Delta x)^2} \left(u_{p}^{n_x-1,j,k,n} - 2 u_{p}^{n_x,j,k,n} 
  + u_{p}^{n_x+1,j,k,n} \right) \right.  \\
  &+\frac{1}{(\Delta y)^2} \left(u_{p}^{n_x,j-1,k,n} - 2 u_{p}^{n_x,j,k,n} + u_{p}^{n_x,j+1,k,n}\right)  \\
  &+\left. \frac{1}{(\Delta z)^2} \left(u_{p}^{n_x,j,k-1,n} - 2 u_{p}^{n_x,j,k,n} + u_{p}^{n_x,j,k+1,n} \right) \right]. 

From the recursions, we get

.. math::

  &a_p \left( u_{p-1}^{n_x-1,j,k,n+1} + u_{p-1}^{n_x,j,k,n+1} - u_{p-1}^{n_x-1,j,k,n} - u_{p-1}^{n_x,j,k,n} \right)    \\
  &+\frac{c \Delta t}{\Delta x} \left(u_{p-1}^{n_x,j,k,n+1} - u_{p-1}^{n_x-1,j,k,n+1} + u_{p-1}^{n_x,j,k,n} 
    - u_{p-1}^  {n_x-1,j,k,n} \right) \nonumber \\
  &+ \sigma_p \frac{\Delta t}{2} \left(u_{p-1}^{n_x-1,j,k,n+1} + u_{p-1}^{n_x,j,k,n+1} + u_{p-1}^{n_x-1,j,k,n} 
    + u_{p-1}^{n_x,j,k,n} \right)  \nonumber \\
  & = \bar{a}_p \left(u_{p}^{n_x-1,j,k,n+1} + u_{p}^{n_x,j,k,n+1} - u_{p}^{n_x-1,j,k,n} - u_{p}^{n_x,j,k,n} \right)  \\
  & -\frac{c \Delta t}{\Delta x} \left(u_{p}^{n_x,j,k,n+1} - u_{p}^{n_x-1,j,k,n+1} + u_{p}^{n_x,j,k,n} 
    - u_{p}^{n_x-1,j,k,n} \right)  \nonumber \\
  & + \bar{\sigma}_p \frac{\Delta t}{2} \left(u_{p}^{n_x-1,j,k,n+1} + u_{p}^{n_x,j,k,n+1} + u_{p}^{n_x-1,j,k,n} 
    + u_{p}^{n_x,j,k,n} \right),  \nonumber 

and

.. math::

  & a_p \left(u_{p-1}^{n_x,j,k,n+1} + u_{p-1}^{n_x+1,j,k,n+1} - u_{p-1}^{n_x,j,k,n} - u_{p-1}^{n_x+1,j,k,n} \right)  \\
  &+\frac{c \Delta t}{\Delta x} \left(u_{p-1}^{n_x+1,j,k,n+1} - u_{p-1}^{n_x,j,k,n+1} + u_{p-1}^{n_x+1,j,k,n} 
    - u_{p-1}^{n_x,j,k,n} \right)  \nonumber \\
  &+ \sigma_p \frac{\Delta t}{2} \left(u_{p-1}^{n_x,j,k,n+1} + u_{p-1}^{n_x+1,j,k,n+1} + u_{p-1}^{n_x,j,k,n} 
    + u_{p-1}^{n_x+1,j,k,n} \right)  \nonumber \\
  &= \bar{a}_p \left(u_{p}^{n_x,j,k,n+1} + u_{p}^{n_x+1,j,k,n+1} - u_{p}^{n_x,j,k,n} - u_{p}^{n_x+1,j,k,n} \right) \\
  & -\frac{c \Delta t}{\Delta x} \left(u_{p}^{n_x+1,j,k,n+1} - u_{p}^{n_x,j,k,n+1} + u_{p}^{n_x+1,j,k,n} 
    - u_{p}^{n_x,j,k,n} \right)  \nonumber \\
  &+ \bar{\sigma}_p \frac{\Delta t}{2} \left(u_{p}^{n_x,j,k,n+1} + u_{p}^{n_x+1,j,k,n+1} + u_{p}^{n_x,j,k,n} 
   + u_{p}^{n_x+1,j,k,n} \right).  \nonumber 

Finally, the termination condition is discretized by

.. math::

  \left(u_{P}^{n_x,j,k,n+1} + u_{P}^{n_x+1,j,k,n+1} - u_{P}^{n_x,j,k,n} - u_{P}^{n_x+1,j,k,n} \right) &  \\
  +\frac{c \Delta t}{\Delta x} \left(u_{P}^{n_x+1,j,k,n+1} - u_{P}^{n_x,j,k,n+1} + u_{P}^{n_x+1,j,k,n} 
  - u_{P}^{n_x,j,k,n} \right) & = 0.

The boundaries on the remaining faces are similar. To deal with edges, we
instead introduce a doubly indexed set of auxiliary variables where the two indexes
correspond to applying the DAB boundary as above in each of the normal directions.
Corners are handled similarly with a triply indexed set of auxiliary variables.

.. rubric:: References

.. bibliography:: zcite.bib
   :list: enumerated
