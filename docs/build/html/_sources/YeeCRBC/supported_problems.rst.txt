******************
Supported Problems
******************

Materials
=========

The library currently requires that materials be isotropic at and beyond the
location of the boundary. It is fine to have arbitrary media in the interior of
the simulation, but any such material should be seperated from the boundary by a 
distance of at least :math:`\delta`.

Mathematically, the reason for this is that inhomogeneities in the materials may
result in waves being (correctly) reflected back into the domain from beyond the
location we place an artificial boundary. The best way to handle these situations,
in general, is to use a larger compuational domain (if possible).

In principle, it is possible to handled layered media. We are considering adding
this ability in a future release.

Initial Conditions
==================

Initial conditions need to be compactly supported (or at least appear to be numerically)
away from the artificial boundaries. That is any initial conditions should be 
zero within :math:`\delta` of the boundary condition.

The issue here is that with non-zero initial conditions, we need a way to 
initialize the auxiliary variables we use in the boundary condition. In general,
we do not know how to initialize the auxiliary variables to arbitrary conditions.

Other Boundary Conditions
=========================

Currently, we support any configuration of rectangular homogeneous Dirichlet, 
homogeneous Neumann, and CRBC/DAB boundaries intersecting at right angles. Note
that it is possible to implement PEC/PMC boundaries under these constraints. We do
not currently support periodic conditions.
