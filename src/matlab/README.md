

<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank">
    <img src="https://github.com/cdslaborg/paramonte/blob/main/img/paramonte.png" alt="ParaMonte: Parallel Monte Carlo and Machine Learning Library" style="width:500px;height:500px;" />
</a>
<br>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=purple&label=release%20version&style=flat-square" alt="GitHub release (latest by date)" /></a>
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a>
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://www.cdslab.org/paramonte/codecov/fortran/2/serial/" target="_blank"><img src="https://img.shields.io/badge/Fortran%20code%20coverage-serial%20:%2077.0%25-brightgreen?style=flat-square" alt="Fortran code coverage - serial" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://www.openhub.net/p/paramonte" target="_blank"><img src="https://img.shields.io/badge/Open%20Hub-stats?color=brightgreen&label=stats&message=Open%20Hub&style=flat-square" alt="stats - Open Hub" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=GitHub%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=PyPI%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://pypistats.org/packages/paramonte" target="_blank"><img src="https://img.shields.io/badge/stats-green?style=flat-square&label=PyPI&labelColor=grey&color=brightgreen" alt="PyPI stats" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/main/src" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/main" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://zenodo.org/record/4076479#.X4Stte17ng4" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4076479.svg" alt="citations and references" /></a>
<a href="https://www.cdslab.org/paramonte/generic/latest/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<a style="border-width:0" href="https://doi.org/10.21105/joss.02741"><img src="https://joss.theoj.org/papers/10.21105/joss.02741/status.svg?style=flat-square" alt="DOI badge" ></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a>
<br><br>
<a href="#paramonte-plain-powerful-parallel-monte-carlo-library">Overview</a> |
<a href="#installation">Installation</a> |
<a href="#dependencies">Dependencies</a> |
<a href="#parallelism">Parallelism</a> |
<a href="#example-usage-instructions">Examples</a> |
<a href="#citing-paramonte">Acknowledgments</a> |
<a href="#license">License</a> |
<a href="#authors-and-contributors">Authors</a>
</div>

ParaMonte: Parallel Monte Carlo and Machine Learning  
====================================================  

ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical density functions of arbitrary dimensions and 
Machine Learning (ML) algorithms for scientific inference, with the design goal of unifying **automation** (of simulations and tasks), 
**user-friendliness** (of algorithms), **accessibility** (from any platform or programming environment), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit https://www.cdslab.org/paramonte  

ParaMonte Design Goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of Monte Carlo and Machine Learning simulations as much as possible to ensure user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing simulation models. 

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, and Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level library implementation to ensure the fastest-possible inferences.  

-   **Parallelizability** of all simulations via shared-memory OpenMP threads and two-sided and one-sided 
    MPI/Coarray distributed communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations with arbitrary digits of precision as requested by the user.  

-   **Comprehensive reporting and post-processing** of simulations and ML tasks, as well as their automatic storage in 
    external files to ensure the simulation results will be understandable and reproducible at any time in the foreseeable future.  


Quick Start
===========


+   To build the library from source or download the prebuilt ParaMonte MATLAB libraries, visit
    [https://www.cdslab.org/paramonte/generic/latest/installation/matlab/](https://www.cdslab.org/paramonte/generic/latest/installation/matlab/).

+   For a quick start with ParaMonte MATLAB examples, visit 
    [this ParaMonte documentation page](https://www.cdslab.org/paramonte/generic/latest/usage/examples/matlab/).

+   For a quick start with all ParaMonte MATLAB functionalities, visit 
    [the ParaMonte MATAB documentation website](https://www.cdslab.org/paramonte/matlab/3/).

The following example code samples a 4-dimensional MultiVariate Normal (MNV) distribution via the ParaDRAM sampler in serial mode,

```matlab
addpath(genpath("./"),"-begin") % change this path to the root directory of paramonte where +pm folder exists.
sampler = pm.sampling.Paradram();
getLogFunc = @(x) -0.5 * sum( x.^2 );
sampler.run ( getLogFunc ...           the objective function
            , 4 ... assume a 4-dimensional objective function
            );
sample = sampler.readSample(); % read the output sample(s).
sample{1}.vis.triplex.lshc3.make() % Make a triplex (corner) plot of the first output sample.
sample{1}.vis.triplex.lshc2.make() % Make a triplex (corner) plot of the first output sample.
```


Installation  
============  


We strongly advise you to follow the guidelines in 
[this documentation page](https://www.cdslab.org/paramonte/generic/1/installation/matlab/) 
to either download the pre-built ParaMonte MATLAB libraries for your platform or to build the library from the source code.  

> Help us tailor the ParaMonte library to your needs through 
> [this single-question survey](https://utaedu.questionpro.com/t/AS6OrZ4xIY).  


Dependencies
============

Aside from an optional MPI runtime library for MPI-parallel simulations, 
the ParaMonte library core has **zero dependency** on external third-party libraries or packages.  



Parallelism
===========

The ParaMonte library relies on three independent parallelism paradigms for parallel applications:  

1.  The Coarray parallelism (only available for select routines in the Fortran programming language).  
2.  The Message Passing Interface (MPI) standard for inter-processor distributed parallelism.  
    +   On **Windows** and **Linux** operating systems, we highly recommend downloading and 
        installing the [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
        which is available to the public free of charge.  
        The [ParaMonte auxil prerelease](https://github.com/cdslaborg/paramonte/releases/tag/auxil)
        exclusively contains the optional Intel MPI runtime libraries for users' convenience.  
        The Intel MPI requires the ParaMonte library to be (or have been prebuilt) using the Intel compilers.  
    +   On **macOS** (where Intel MPI is unavailable) and **Linux** operating systems, users can install
        either [Open-MPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/) MPI runtime 
        libraries depending on the prebuilt version of the ParaMonte library that you have 
        downloaded or the configuration with which you intend to build the library.  

**Optional Dependency Compatibility**  

| Dependency                        | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| Intel MPI (IMPI) >= 2021.8        | ✅              | ✅            | ✅             | ❌            |  
| MPICH MPI (MMPI) >= 3             | ❌              | ✅            | ✅             | ✅            |  
| OpenMPI (OMPI) >= 4               | ❌              | ✅            | ✅             | ✅            |  

Read the [`install.md` installation instruction notes](install.md) to 
learn more about the library's optional parallel build configurations.<br>

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).


+   The serial version of the ParaMonte MATLAB library has **NO external library dependencies**.  
+   The MPI-parallel version requires a ParaMonte-compatible MPI runtime library.  
+   The new thread-based parallelism features of the ParaMonte MATLAB library require the MATLAB Parallel Toolbox.  
+   The MATLAB Parallel Toolbox is usually automatically shipped with MATLAB software.  


Example Usage Instructions
==========================

+   For complete clear organized, up-to-date instructions on the build process and the 
    installation of the ParaMonte library, visit the [ParaMonte MATLAB examples guidelines](https://www.cdslab.org/paramonte/generic/latest/usage/examples/matlab/).  


Citing ParaMonte
================

The ParaMonte library is an honor-ware and its currency is acknowledgment and citations. 

If you use ParaMonte or any ideas from this software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/generic/latest/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work)
to access the PDF version of these files free of charge.  


License
=======

The majority of the ParaMonte library routines are distributed under the permissive 
[MIT License](LICENSE.md) plus acknowledgments that are detailed in [LICENSE.md](./LICENSE.md).  

However, there are occasionally modernized and extended routines from external projects, particularly within the ParaMonte Fortran library. 
These routines may have a license that does not fully overlap with the default ParaMonte library license [LICENSE.md](./LICENSE.md). 
In such cases, the documentation and source code of the specific routine explicitly list the original license information. 

This is free software, so help us keep it freely available to the public by acknowledging the usage and contributing to it. 
If you have questions or concerns about the license, please get in touch with the project's lead (shahmoradi@utexas.edu). 


Authors and Contributors
========================

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)
    +   Astrophysicist/bioinformatician by training,
    +   Ph.D. in computational physics/bioinformatics, The University of Texas at Austin,
    +   Currently a faculty member of the Data Science and Physics programs at UT Arlington,
    +   Contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")

+   [Joshua Osborne](https://www.cdslab.org/people/#joshua-alexander-osborne)
    +   Physicist / computational data scientist by training,
    +   Currently a Physics PhD candidate at UT Arlington,
    +   Contact: [joshuaalexanderosborne@gmail.com](mailto:"joshuaalexanderosborne@gmail.com")

+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)
    +   PhD in Astronomy,
    +   PhD in Space Physics,
    +   Deep philosophical thinker,
    +   Currently at NASA Goddard Space Flight Center,
    +   Contact: [bagheri.fateme@gmail.com](mailto:"bagheri.fateme@gmail.com")


Contributing to ParaMonte
=========================

Please read the generic contribution guidelines listed in [CONTRIBUTING.md](./CONTRIBUTING.md).


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")
