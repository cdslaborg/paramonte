[![ParaMonte - Plain Powerful Parallel Monte Carlo Library](https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png)](https://www.cdslab.org/paramonte/)  
  
<a href="https://github.com/cdslaborg/paramonte/blob/master/LICENSE" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a> 
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a> 
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a> 
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/master/src/interface" target="_blank"><img src="https://img.shields.io/badge/interface-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=pypi%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/master/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/master" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a> 
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), **user-friendliness** (of the library), **accessibility** (from multiple programming environments), **high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  

**For a quick start with some MATLAB Live Script examples, visit [this ParaMonte documentation page](https://www.cdslab.org/paramonte/notes/examples/matlab/mlx/). The corresponding example source files (the `*.mlx` files) can be downloaded from the [paramontex GitHub repository](https://github.com/cdslaborg/paramontex/tree/master/MATLAB/mlx), a repository dedicated to the ParaMonte library examples.**  

ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

+   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  
+   **Interoperability** of the core library with as many programming languages as currently possible, including C/C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  
+   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  
+   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray communications while requiring zero-parallel-coding efforts by the user.  
+   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  
+   **Fully-deterministic reproducibility** and automatically-enabled restart functionality for all simulations up to 16 digits of precision as requested by the user.  
+   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

Installation  
============  

The latest release of the ParaMonte MATLAB library can be downloaded from the release page of the library's repository on GitHub:  

[https://github.com/cdslaborg/paramonte](https://github.com/cdslaborg/paramonte)  

Alternatively, you can build the library from the source in the GitHub repository of the project ([https://github.com/cdslaborg/paramonte](https://github.com/cdslaborg/paramonte)). For instructions, please visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

Dependencies  
============  

The serial version of the MATLAB interface to the ParaMonte kernel library has no external library dependencies.  

Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). These automatic checks and installations happen when you download, and use the library, for the first time, in your MATLAB environment. If the automatic installation is unsuccessful, you can also install the libraries manually on your system. On Windows and Linux operating systems, we highly recommend downloading and installing the [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), which is available to the public free of charge. On macOS, we recommend [Open-MPI](https://www.open-mpi.org/) since the Intel MPI library does not support macOS. For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

+   On Windows and Linux operating systems, we highly recommend downloading and installing the [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), which is available to the public free of charge.  
+   On macOS, we recommend [Open-MPI](https://www.open-mpi.org/) since the Intel MPI library does not support macOS. For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

Citing ParaMonte  
================  

As per the license terms, we kindly ask you to cite the ParaMonte library if you use the library or any parts of it in your research, education, or software development. Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) for the full citation and reference information.  

Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi):  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

- [Shashank Kumbhare](https://www.cdslab.org/people/#shashank-kumbhare):  
    - physicist / Computational Data Scientist by training,  
    - currently a UTA Physics member,  
    - contact: [shashankkumbhare8@gmail.com](mailto:"shashankkumbhare8@gmail.com")  

License  
=======  
  
[GNU Lesser General Public License](https://github.com/cdslaborg/paramonte/blob/master/LICENSE)  
  
Example usage instructions  
==========================  

+   **Install a MATLAB >2017b distribution**, preferably, the the latest MATLAB. Note that ParaMonte MATLAB library have been tested only with MATLAB version 2018b and newer.  

+   **Optionally install a compatible MPI library** (or let the ParaMonte library take care of the installation when you call the library for the first time). For parallel simulations (via MPI), you will need an MPI library already installed on your system. If you choose to install the library by yourself, we recommend the Intel MPI library which is available for free from the Intel website. On macOS, the OpenMPI library can be used in place of the Intel MPI library which currently does not support macOS.  

+   **Calling the ParaMonte library for the first time**  
    +   **Windows**
        Nothing special needs to be done. You are all set! Follow the instructions below on how to run your ParaMonte-enabled simulations.  
    +   **Linux/macOS**
        Nothing special needs to be done. You are all set! Follow the instructions below on how to run ParaMonte-enabled simulations.  

+   **Running the ParaMonte simulations**  
    +   Open you.  
    +   Suppose your mathematical objective function is a multivariate Normal distribution as implemented in this [logfunc.py](https://raw.githubusercontent.com/cdslaborg/paramonte/master/example/mvn/Python/logfunc.py) file.  
    +   For **serial** simulations, download this example generic serial [main.py](https://raw.githubusercontent.com/cdslaborg/paramonte/master/example/main.py) Python main file and save it in the same folder containing the `logfunc.py` file that you downloaded in the above. Then, simply type the name of the Python main script, `python main.py` on the Bash terminal or the Anaconda command line.  
    +   For **parallel** simulations, download this example generic serial [main_mpi.py](https://raw.githubusercontent.com/cdslaborg/paramonte/master/example/main_mpi.py) Python main file and save it in the same folder containing the `logfunc.py` file that you downloaded in the above. Then, simply invoke the MPI launcher followed by the name of the Python main script on the Bash terminal, similar to the following,  
        +   on Windows (within the Anaconda command line or a terminal that recognizes both `mpiexec` and `python` software),  
            ```  
            mpiexec -localonly -n 3 python main_mpi.py
            ```  
        +   on macOS or Linux (within a Bash terminal),  
            ```  
            mpiexec -n 3 python main_mpi.py
            ```  
        Here, the parallel simulations are performed on 3 processes. Change the number 3 to any number of processes you wish to use, but do not go beyond the maximum number of physical processes available on your system, otherwise, it will only degrade the performance of your parallel simulations. For example, if you are running the parallel simulation on a personal quad-cores laptop, set the number of processes to either 3 or 4 at most.  
    +   Enjoy the unification of simplicity, efficiency, and parallelism in Monte Carlo simulations!  
    +   The ParaMonte library samplers are extremely versatile with many adjustable input parameters. To learn about the advanced features of the ParaMonte routines, visit: https://www.cdslab.org/paramonte  
