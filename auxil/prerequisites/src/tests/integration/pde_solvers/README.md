[This document is formatted with GitHub-Flavored Markdown.   ]:#
[For better viewing, including hyperlinks, read it online at ]:#
[https://github.com/sourceryinstitute/OpenCoarrays/edit/master/src/tests/integration/pde_solvers/README.txt]:#

Partial Differential Equation (PDE) Solvers
===========================================

This directory contains three PDE solvers listed here in order from simplest to most complex:

* A one-dimensional (1D) finite-difference, unsteady [heat equation solver],
* A 1D finite-difference, unsteady, nonlinear [Burgers equation solver], and
* A three-dimensional (3D), unsteady, spectral [Navier-Stokes equation solver].

The first two solvers contain correctness checks that result in the printing of the 
message "Test passed" before terminating a correct execution.  For more details on the 
heat equation solver please view the [Sourcery Institute] [tutorial videos] online.  

For more details on the Burgers solver, please see Chapter 12 of the textbook 
[Scientific Sofware Design] or the open-access journal article 
"[High Performance Design Patterns for Modern Fortran]."  The [coarrayBurgers] 
subdirectory includes a [run.sh] launch script that works inside the open-source 
Linux virtual machine available in the Sourcery Institute [store].

The launch script instruments the Burgers solver for performance analysis using the 
open-source Tuning and Analysis Utilities ([TAU]) package.  The instrumented Burgers 
solver has been demonstrated to execute with 87% parallel efficiency on 16,384
cores in weak scaling when compiled with the Cray Compiler Environment.  For
new scalabiliby studies, it is important to run problems of sufficient size.  For
performance and complexity comparisons, an MPI version of the Burgers solver is
in the [performance] directory at the same level as the current directory.

The Navier-Stokes solver uses Fourier-spectral methods and Runge-Kutta time advancement
to simulate the evolution of statistically homogeneous turbulent flow in a 3D box with 
periodic boundary conditions.  For performance and complexity comparisons, the 
[navier-stokes] subdirectory contains both a  Message Passing Interface (MPI) version 
and a coarray Fortran (CAF) version of the same solution algorithm.

[heat equation solver]: ./coarrayHeatSimplified
[Burgers equation solver]: ./coarrayBurgers
[Navier-Stokes equation solver]: ./navier-stokes
[Sourcery Institute]: http://www.sourceryinstitute.org
[tutorial videos]: http://www.sourceryinstitute.org/videos
[Scientific Sofwtware Design]: http://www.cambridge.org/rouson
[High Performance Design Patterns for Modern Fortran]: http://www.hindawi.com/journals/sp/2015/942059/
[store]: http://www.sourceryinstitute.org/store
[coarrayBurgers]: ./coarrayBurgers
[run.sh]: ./coarrayBurgers/run.sh
[TAU]: http://tau.uoregon.edu
[navier-stokes]: ./navier-stokes
[performance]: ../../performance
