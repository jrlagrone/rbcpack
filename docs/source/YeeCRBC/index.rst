.. For component pages, the structure should be the following:
   1) The component name should be at the top, underlined with "==="
      This will ensure that it ends up in the "components" drop down menu
      in the top bar.

   2) Subsequent headings should be used in the following order of hierarchy:
      ---
        ^^^
          """
      Any heading in this file underlined with "===", "---", "^^^" will be included
      in the side bar navigation nested based on the above. "===" should only
      be used for component titles for consistancy 

    3) This file needs to be linked to the rest of the project by including it in
       the toctree list in the top-level index.rst file

Yee/FDTD
========

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

Downloads and Respository
-------------------------

The initial release download at `yee_crbc_lib_1.0.0.tar.gz <yee_crbc_lib_1.0.0.tar.gz>`_

The code is also hosted in a read only repository that may contain bug fixes and
improvements between releases. `RBCPack Bitbucket Reposistory <https://bitbucket.org/rbcpack/rbcpack>`_

Licensing
^^^^^^^^^

The code is released under the LGPL v3. See http://www.gnu.org/licenses/lgpl-3.0.en.html.
If a different license is required, please contact us.

Documentation
-------------
User Guides
^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   install
   supported_problems
   gen_use
   examples

In addition, Doxygen generated documentation is available: 
`Doxygen Documentation <./doxygen/html/index.html>`_

Known Issues
^^^^^^^^^^^^

There are currently no known issues.

Theory, Results, and Publications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   theory_overview
   matlab_examples
   results
   pubs

Acknowledgments
---------------

This work was supported in part by ARO Grant W911NF-09-1-0344 and DoD STTR A11A-015-0067. 
All conclusions are those of the authors and do not necessarily reflect
the views of the ARO or DoD.
