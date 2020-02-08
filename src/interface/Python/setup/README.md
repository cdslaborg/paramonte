[![ParaMonte - Plain Powerful Parallel Monte Carlo Library](https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png)](https://cdslab.org/paramonte/)  
  
<br>
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), **user-friendliness** (of the library), **accessibility** (from multiple programming environments), **high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  

ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-  **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest-level of user-friendliness of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  
-  **Interoperability** of the core library with as many programming languages as currently possible, including C/C++, Fortran, Python, with ongoing efforts to support other popular programming languages.  
-  **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  
-  **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray communications while requiring zero-parallel-coding efforts from the user.  
-  **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  
-  **Fully-deterministic reproducibility** and automatically-enabled restart functionality for all simulations up to 16 digits of precision as requested by the user.  
-  **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

Installation  
============  

The latest release of ParaMonte can be installed from PyPI using `pip`:  

    pip install --user paramonte  

Alternatively, you can build the library from the source in the GitHub repository of the project ([https://github.com/cdslaborg/paramonte](https://github.com/cdslaborg/paramonte)). For instructions, please visit: [cdslab.org/pm](https://github.com/cdslaborg/paramonte)  

Dependencies  
============  

The Python interface of ParaMonte depends on several Python standard modules, as well as third-party libraries. These include numpy, scipy, pandas, matplotlib, and seaborn. The last two (plotting) libraries are only used for the post-processing of simulation results and are therefore not needed if you do not plan to use the post-processing functionalities of ParaMonte. If you have a recent version of Anaconda Python distribution installed on your system, then all of the dependencies already exist and are automatically installed on your system.  

Parallelism  
===========  

ParaMonte relies on Message Passing Interface (MPI) standard for inter-processor communications. To run a parallel simulation you will have a compatible MPI library installed on your system. In most cases, ParaMonte can automatically install the required MPI libraries on your system (with your permission), when you download, install, and `import paramonte` into your Python environment for the first time. If the automatic installation is unsuccessful, you can also install the libraries manually on your system. On Windows and Linux operating systems, we highly recommend downloading and installing [Intel MPI library](https://software.intel.com/en-us/mpi-library), which is available free of charge. On macOS, we recommend [MPICH](https://www.mpich.org/downloads/) since the Intel MPI library does not support macOS.  

Citing ParaMonte  
================  

The paper for this package is already in publication process. Until then, if you use this package, we kindly ask you to cite:  

-  Amir Shahmoradi (2013). A Multivariate Fit Luminosity Function and World Model for Long Gamma-Ray Bursts. The Astophysical Journal (ApJ) 766:111-133 ([PDF link](https://www.cdslab.org/pubs/Shahmoradi_2013a.pdf)).  

Authors  
=======  

- [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi):  
    - astrophysicist / bioinformatician by training (and a science-lover in general),  
    - Ph.D. in computational physics / bioinformatics from the University of Texas at Austin,  
    - currently a faculty member of Physics and Data Science at The University of Texas at Arligton,  
    - with teaching / research experience/background in computational and data sciences, statistics, data analysis and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, protein dynamics and interactions),  
    - contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

- [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri):  
    - physicist / cosmologist by training,  
    - currently a UTA Physics member,  
    - deep philosophical thinker,  
    - contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  

License  
=======  

[GNU Lesser General Public License](https://github.com/cdslaborg/paramonte/blob/master/LICENSE.txt)
