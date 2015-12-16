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
such as Maxwell's equations or accoustics.

These routines are available in the `RBCPack Bitbucket Reposistory <https://bitbucket.org/rbcpack/rbcpack>`_
or can be downloaded directly `optimal_cosines.tar.gz <optimal_cosines.tar.gz>`_

