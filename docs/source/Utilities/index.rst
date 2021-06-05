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

Utilities
=========

The following are some useful routines that may be useful for writing a custom
implementation of CRBC/DAB boundaries

Optimal Cosines
---------------

These files compute optimal parameters for CRBC/DAB boundaries for the scalar
wave equation. It is appropriate to use these values for most homogeneous problems
such as Maxwell's equations or accoustics. There are two functions provided:
*optimal_cosines* and *optimal_cosinesP*. *optimal_cosines* computes optimal
cosine parameters and the optimal number of recursions based on an user provided
error tolerance. *optimal_cosines* computes the optimal cosine parameters based
on an user provided number of recursions. The usage is documented in the header
file optimal_cosines.h.

These routines are available in the `RBCPack github Reposistory <https://github.com/jrlagrone/rbcpack>`_
or can be downloaded directly :download:`optimal_cosines.tar.gz <optimal_cosines.tar.gz>`.

