Welcome to the Radiation Boundary Condition Package
===================================================

The purpose of this library is to provide accurate and efficient radiation and
absorbing boundary conditions to truncate the computational domain for the
simulation of a variety of wave propagation problems in the time domain. 
The plan is for RBCPack to contain stand-alone components that implement complete
radition boundary conditions and/or double absorbing boundary layers for a variety
of popular time-domain wave propagation solvers.


Components
----------

Currently, the following components are available:

.. We do the toctree to generate the list of pages that will appear in the 
   "component" drop down menu in the top navigation bar. Subsequent components
   should be added to this list to ensure that sphinx can find them ...

.. toctree::
   :hidden:
   :maxdepth: 1

   YeeCRBC/index
   Utilities/index

.. Now we make some "fancy" boxes for each of the components

.. raw:: html

  <div class="list-group">
    <a href="YeeCRBC/index.html" class="list-group-item">
      <h4 class="list-group-item-heading">Yee/FDTD</h4>
      <p class="list-group-item-text">
          This component provides complete radiation boundary conditions in the
          for of a double absorbing boundary for Yee/FDTD solvers. </p>
    </a>
    <a href="Utilities/index.html" class="list-group-item">
      <h4 class="list-group-item-heading">Utilities</h4>
      <p class="list-group-item-text">
         A collection of useful routines that may aid in the development of other
         boundary implementations.</p>
    </a>
  </div>

Other components are under development and will be released as they are completed.
Components to handle (un)structured discontinuous Galerkin methods and methods
based on modified equations are nearing completion and should be released in the
near future. Components for Hermite methods are in the early stages of development.

Downloading
-----------

The code is also hosted in a read only repository: 
`RBCPack Bitbucket Reposistory <https://bitbucket.org/rbcpack/rbcpack>`_.
Alternatively, the current release of each component can be downloaded from their
respective pages.

The code is released under the LGPL v3. See http://www.gnu.org/licenses/lgpl-3.0.en.html.
If a different license is required, please contact us.



Acknowledgments
^^^^^^^^^^^^^^^

This work was supported in part by ARO Grant W911NF-09-1-0344 and DoD STTR A11A-015-0067. 
All conclusions are those of the authors and do not necessarily reflect
the views of the ARO or DoD.

