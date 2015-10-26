.. Yee CRBC Library documentation master file, created by
   sphinx-quickstart on Tue Jan 13 12:32:34 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Yee CRBC Library
===========================

The purpose of this library is to provide accurate and efficient radiation and
absorbing boundary conditions to truncate the computational domain for the
simulation of wave propagation problems in the time domain. This particular component
of the library is specialized to second order finite difference time domain (FDTD) 
simulations, in particular, the Yee scheme.

This library utilizes a local boundary condition formulation called Complete 
Radiation Boundary Conditions (CRBCs). The CRBCs have been implemented in a
thin layer to deal with the staggered space-time grid common in FDTD discretizations
and we refer to this construction as a Double Absorbing Boundary (DAB) layer.
The major benefit of the CRBC/DAB boundary construction is that there exists a 
clear notion of convergence and it is possible to choose optimal parameters `a priori`.

Acknowledgments:
----------------

This work was supported in part by ARO Grant W911NF-09-1-0344 and DoD STTR A11A-015-0067. 
All conclusions are those of the authors and do not necessarily reflect
the views of the ARO or DoD.

Software
========

Licensing
---------

We have released the code under the LGPL v3. See http://www.gnu.org/licenses/lgpl-3.0.en.html.

Download
--------

The initial release candidate can be download at `yee_crbc_lib_1.0.0_rc.tar.gz <yee_crbc_lib_1.0.0_rc.tar.gz>`_

Repository
----------

A read only repository will be available in the near future to provide updates.

Documentation
=============

Known Issues
------------

There are currently no known issues.

User Guides
-----------

.. toctree::
   :maxdepth: 2

   install
   gen_use
   examples

In addition to Doxygen generate documentation is available online <here>

Theory, Results and Publications
================================

.. toctree::
   :maxdepth: 2

   theory_overview
   results
   pubs
