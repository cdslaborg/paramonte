[![ParaMonte: Plain Powerful Parallel Monte Carlo Library](https://www.cdslab.org/paramonte/images/paramonte.png)](https://www.cdslab.org/paramonte)  
  
[![GitHub](https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square)](https://github.com/cdslaborg/paramonte/blob/master/LICENSE) 
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square)](https://github.com/cdslaborg/paramonte/releases) 
![GitHub Release Date](https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square) 
[![PyPI - release version](https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square)](https://pypi.org/project/paramonte/) 
![PyPI - Status](https://img.shields.io/pypi/status/paramonte?style=flat-square) 
[![supported languages](https://img.shields.io/badge/interface-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square)](https://github.com/cdslaborg/paramonte/tree/master/src/interface) 
![GitHub All Releases](https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square) 
[![PyPI - Downloads](https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=pypi%20downloads&style=flat-square)](https://pypi.org/project/paramonte/) 
![GitHub repo size](https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square) 
![GitHub language count](https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square) 
[![GitHub commit activity](https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square)](https://github.com/cdslaborg/paramonte/graphs/contributors) 
[![citations and references](https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square)](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work)  

[![Twitter](https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte)](https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte)
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), **user-friendliness** (of the library), **accessibility** (from multiple programming environments), **high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-  **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  
-  **Interoperability** of the core library with as many programming languages as currently possible, including C/C++, Fortran, Python, with ongoing efforts to support other popular programming languages.  
-  **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  
-  **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray communications while requiring zero-parallel-coding efforts by the user.  
-  **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  
-  **Fully-deterministic reproducibility** and automatically-enabled restart functionality for all simulations up to 16 digits of precision as requested by the user.  
-  **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

Installation  
============  

The ParaMonte library installation/build process is fully automated for all of the supported programming languages. The pre-built libraries are also available on [the release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases). Each prebuilt ParaMonte library automatically ships with a full-fledged set of example codes and build scripts.  

For more information and quick-start in the programming language of your choice, visit the [ParaMonte library homepage](https://www.cdslab.org/paramonte).  

Dependencies  
============  

Beyond an optional MPI runtime library for parallel simulations, the ParaMonte kernel has **zero dependency** on external third-party libraries or packages.  

Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). These automatic checks and installations happen when you download and build the library for the first time on your system. If the automatic installation is unsuccessful, you can also install the libraries manually on your system.  
+   On Windows and Linux operating systems, we highly recommend downloading and installing the [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), which is available to the public free of charge.  
+   On macOS, we recommend [Open-MPI](https://www.open-mpi.org/) since the Intel MPI library does not support macOS. For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

Citing ParaMonte  
================  

As per the license terms, we kindly ask you to cite the ParaMonte library if you use the library or any parts of it in your research, education, or software development. Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) for the full citation and reference information.  

Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)  
    +   physicist / cosmologist by training,  
    +   currently a UTA Physics member,  
    +   deep philosophical thinker,  
    +   contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  

+   [Shashank Kumbhare](https://www.cdslab.org/people/#shashank-kumbhare)  
    +   physicist / Computational Data Scientist,  
    +   currently a UTA Physics member,  
    +   contact: [shashankkumbhare8@gmail.com](mailto:"shashankkumbhare8@gmail.com")  

+   [Joshua Osborne](https://www.cdslab.org/people/#joshua-alexander-osborne)  
    +   physicist / Computational Data Scientist by training,  
    +   currently a UTA Physics member,  
    +   contact: [joshuaalexanderosborne@gmail.com](mailto:"joshuaalexanderosborne@gmail.com")  

License  
=======  
  
[GNU Lesser General Public License](https://github.com/cdslaborg/paramonte/blob/master/LICENSE)  
  
Example usage instructions  
==========================  

+   For complete organized up-to-date instructions, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

+   For a quick look into *language-specific* README.md instructions, visit:  
    +   **C/C++**: [https://github.com/cdslaborg/paramonte/tree/master/src/interface/C](https://github.com/cdslaborg/paramonte/tree/master/src/interface/C)  
    +   **Fortran**: [https://github.com/cdslaborg/paramonte/tree/master/src/interface/Fortran](https://github.com/cdslaborg/paramonte/tree/master/src/interface/Fortran)  
    +   **MATLAB**: [https://github.com/cdslaborg/paramonte/tree/master/src/interface/MATLAB](https://github.com/cdslaborg/paramonte/tree/master/src/interface/MATLAB)  
    +   **Python**: [https://github.com/cdslaborg/paramonte/tree/master/src/interface/Python](https://github.com/cdslaborg/paramonte/tree/master/src/interface/Python)  

