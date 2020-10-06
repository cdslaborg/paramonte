
  
<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank"><img src="https://raw.githubusercontent.com/shahmoradi/paramonte/gh-pages/images/paramonte.png" alt="ParaMonte: Plain Powerful Parallel Monte Carlo Library" /></a>
<br><br>
<a href="https://github.com/cdslaborg/paramonte/blob/master/LICENSE.md" target="_blank"><img src="https://img.shields.io/github/license/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub" /></a>  
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/v/release/cdslaborg/paramonte?color=orange&label=kernel%20release&style=flat-square" alt="GitHub release (latest by date)" /></a> 
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/github/release-date/cdslaborg/paramonte?color=orange&style=flat-square" alt="GitHub Release Date" /></a> 
<a href="https://pypi.org/project/paramonte/" target="_blank"><img src="https://img.shields.io/pypi/v/paramonte?color=orange&label=pypi%20release&style=flat-square" alt="PyPI - release version" /></a> 
<a href="https://github.com/cdslaborg/paramonte/releases" target="_blank"><img src="https://img.shields.io/pypi/status/paramonte?style=flat-square" alt="PyPI - Status" /></a>
<a href="https://lgtm.com/projects/g/cdslaborg/paramonte/?mode=list" target="_blank"><img src="https://img.shields.io/lgtm/grade/python/github/cdslaborg/paramonte?label=code%20quality&style=flat-square&color=brightgreen" alt="LGTM Grade" /></a>
<a href="https://github.com/cdslaborg/paramonte/issues" target="_blank"><img src="https://img.shields.io/github/issues/cdslaborg/paramonte?style=flat-square" alt="GitHub issues" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/master/src/interface" target="_blank"><img src="https://img.shields.io/badge/available%20in-C%20%2F%20C%2B%2B%20%2F%20Fortran%20%2F%20MATLAB%20%2F%20Python-brightgreen?style=flat-square" alt="supported languages" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/traffic" target="_blank"><img src="https://img.shields.io/github/downloads/cdslaborg/paramonte/total?color=brightgreen&label=kernel%20downloads&style=flat-square" alt="GitHub All Releases" /></a>
<a href="https://libraries.io/pypi/paramonte" target="_blank"><img src="https://img.shields.io/pypi/dm/paramonte?color=brightgreen&label=pypi%20downloads&style=flat-square" alt="PyPI - Downloads" /></a>
<a href="https://www.mathworks.com/matlabcentral/fileexchange/78946-paramonte" target="_blank"><img src="https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg" alt="View ParaMonte on File Exchange" /></a>  
<a href="https://github.com/cdslaborg/paramonte/" target="_blank"><img src="https://img.shields.io/github/repo-size/cdslaborg/paramonte?style=flat-square" alt="GitHub repo size" /></a>
<a href="https://github.com/cdslaborg/paramonte/tree/master/src/interface" target="_blank"><img src="https://img.shields.io/github/languages/count/cdslaborg/paramonte?style=flat-square" alt="GitHub language count" /></a>
<a href="https://github.com/cdslaborg/paramonte/graphs/contributors" target="_blank"><img src="https://img.shields.io/github/commit-activity/y/cdslaborg/paramonte?style=flat-square" alt="GitHub commit activity" /></a>
<a href="https://github.com/cdslaborg/paramonte/commits/master" target="_blank"><img src="https://img.shields.io/github/last-commit/cdslaborg/paramonte?color=blue&style=flat-square" alt="GitHub last commit" /></a>
<a href="https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work" target="_blank"><img src="https://img.shields.io/badge/reference-%20%09arXiv%3A1209.4647-blueviolet?style=flat-square" alt="citations and references" /></a>
<a href="https://ascl.net/2008.016" target="_blank"><img src="https://img.shields.io/badge/ascl-2008.016-blue.svg?colorB=262255" alt="ascl:2008.016" /></a>
<br><br>
<a href="https://twitter.com/intent/tweet?text=ParaMonte%20-%20Plain%20Powerfull%20Parallel%20Monte%20Carlo%20Library:&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" target="_blank"><img src="https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2Fcdslaborg%2Fparamonte" alt="Twitter" /></a> 
</div>
  
  
ParaMonte: Plain Powerful Parallel Monte Carlo Library
======================================================
  
ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical objective functions 
of arbitrary-dimensions, in particular, the posterior distributions of Bayesian models in data science, Machine Learning, 
and scientific inference, with the design goal of unifying the **automation** (of Monte Carlo simulations), 
**user-friendliness** (of the library), **accessibility** (from multiple programming environments), 
**high-performance** (at runtime), and **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: https://www.cdslab.org/paramonte  
  
  
ParaMonte design goals  
======================  

ParaMonte has been developed while bearing the following design goals in mind:  

-   **Full automation** of all Monte Carlo simulations to the highest levels possible to ensure the highest level of user-friendliness 
    of the library and minimal time investment requirements for building, running, and post-processing of simulation models.  

-   **Interoperability** of the core library with as many programming languages as currently possible, 
    including C, C++, Fortran, MATLAB, Python, with ongoing efforts to support other popular programming languages.  

-   **High-Performance** meticulously-low-level implementation of the library to ensure the fastest-possible Monte Carlo simulations.  

-   **Parallelizability** of all simulations via two-sided and one-sided MPI/Coarray 
    communications while requiring zero-parallel-coding efforts by the user.  

-   **Zero-dependence** on external libraries to ensure hassle-free ParaMonte library builds and ParaMonte simulation runs.  

-   **Fully-deterministic reproducibility** and automatically-enabled restart functionality 
    for all simulations up to 16 digits of precision as requested by the user.  

-   **Comprehensive-reporting and post-processing** of each simulation and its results, as well as their automatic storage in 
    external files to ensure the simulation results will be comprehensible and reproducible at any time in the distant future.  

  
Quick start  
===========  

For a quick start with some MATLAB Live Script examples, visit [this ParaMonte documentation page](https://www.cdslab.org/paramonte/notes/examples/matlab/mlx/). 
The corresponding example source files (the `*.mlx` files) can be downloaded from the [paramontex GitHub repository](https://github.com/cdslaborg/paramontex/tree/master/MATLAB/mlx), 
a repository dedicated to the ParaMonte library examples.  

The following example code samples a 4-dimensional MultiVariate Normal (MNV) distribution via the ParaDRAM sampler in serial mode,  

```matlab  
addpath(genpath("./"),"-begin") % change this path to the root directory of paramonte
pm = paramonte();
pmpd = pm.ParaDRAM();
getLogFunc = @(x) -0.5 * sum( x.^2 );
pmpd.runSampler ( 4 ... assume a 4-dimensional objective function
                , getLogFunc ...           the objective function
                );
```  

To learn about the post-processing and visualization tools of the `ParaMonte::MATLAB` library, visit [this this documentation page](https://www.cdslab.org/paramonte/notes/examples/matlab/mlx/).  

  
Installation  
============  

The latest release of the ParaMonte MATLAB library can be downloaded from the release page of the library's repository on GitHub:  

[https://github.com/cdslaborg/paramonte/releases/latest/](https://github.com/cdslaborg/paramonte/releases/latest/)  

Alternatively, you can build the library from the source in the GitHub repository of the project: https://github.com/cdslaborg/paramonte  
For instructions, please visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

  
Dependencies  
============  

The serial version of the ParaMonte MATLAB library kernel library has **NO external library dependencies**. 
The parallel version requires an MPI runtime library.  

  
Parallelism  
===========  

The ParaMonte library relies on the Message Passing Interface (MPI) standard for inter-processor communications. 
To run a parallel simulation, you will have to have a compatible MPI runtime library installed on your system. 
In most cases, ParaMonte will automatically install the required missing libraries on your system (with your permission). 
These automatic checks and installations happen when you download and install or use the library on your system, for the first time. 
If the automatic installation is unsuccessful, you can also install the libraries manually on your system:  

+   On **Windows** and **Linux** operating systems, we highly recommend downloading and installing the 
    [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), 
    which is available to the public free of charge.  
+   On **macOS**, we recommend [Open-MPI](https://www.open-mpi.org/) since the Intel MPI library does not support macOS. 

For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

  
Citing ParaMonte  
================  

The ParaMonte library is an honor-ware and its currency is acknowledgment and citations.  
  
As per the ParaMonte library license agreement terms, if you use any parts of this library for any purposes, 
kindly acknowledge the use of the ParaMonte library in your work (education/research/industry/development/...) 
by citing the ParaMonte library's main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) 
to access the PDF version of these files free of charge.  

  
License  
=======  

[MIT License](https://github.com/cdslaborg/paramonte/blob/master/LICENSE.md)  

**What does this license mean?**  

Essentially, all we are asking from the users or developers is to  

>   explicitly acknowledge the use of this library or any concepts or parts of it in their education, research, or software (free or commercial).  

This is a free software, so help us keep it freely available to the public by redistributing the library and contributing to it. 
If you have questions or concerns about the license, do not hesitate to contact us (shahmoradi@utexas.edu).  

  
Authors and contributors  
========================  

+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    +   astrophysicist/bioinformatician by training (and a science-lover in general),  
    +   Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    +   currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    +   with teaching/research experience/background in computational and data sciences, statistics, 
        data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, 
        high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, 
        biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, 
        protein dynamics, and interactions),  
    +   contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+   [Shashank Kumbhare](https://www.cdslab.org/people/#shashank-kumbhare)  
    +   physicist / Computational Data Scientist,  
    +   currently a UTA Physics member,  
    +   contact: [shashankkumbhare8@gmail.com](mailto:"shashankkumbhare8@gmail.com")  

  
Example usage instructions  
==========================  

+   **Install a MATLAB >2017b distribution**, preferably, the the latest MATLAB. 
    Note that ParaMonte MATLAB library have been tested only with MATLAB version 2018b and newer.  

+   **Optionally install a compatible MPI library** (or let the ParaMonte library take care of the installation 
    when you call the library for the first time). For parallel simulations (via MPI), you will need an MPI library 
    already installed on your system. If you choose to install the library by yourself, we recommend the Intel MPI 
    library which is available for free from the Intel website. On macOS, the OpenMPI library can be used 
    in place of the Intel MPI library which currently does not support macOS.  

+   **Calling the ParaMonte library for the first time**  
    +   **Windows**  
        Nothing special needs to be done. You are all set! Follow the instructions below on how to run your ParaMonte-enabled simulations.  
    +   **Linux/macOS**  
        For the first time, before running any ParaMonte-enabled simulations, we **highly recommend** you to run MATLAB as an 
        administrator on your system if you can (for example, if you own the system). This can be done from a Bash terminal via,  
        ```bash  
        sudo matlab
        ```  
        If `matlab` is not recognized on your command line as an application, seek help from 
        [this ParaMonte documentation page](https://www.cdslab.org/paramonte/notes/troubleshooting/bash-matlab-command-not-found/). 
        Once the MATLAB interactive environment opens, navigate to the root folder of the ParaMonte library (where the LICENSE file exists) 
        and call the ParaMonte library for the first time via the following commands (simply type the commands on the MATLAB command prompt),  
        ```matlab  
        addpath(genpath("./"),"-begin"); % add the ParaMonte library directories to MATLAB's list of search paths.
        pm = paramonte(); % instantiate an object of class paramonte.
        pm.verify(); % verify the integrity of the ParaMonte library on your system.
        ```  
        If needed, follow any extra instructions provided by the library on your MATLAB command prompt. 
        Before beginning to use the ParaMonte library, we strongly recommend that you close your current MATLAB 
        session and the Bash terminal from which you initiated the MATLAB session in **sudo** mode, **entirely**. 
        Then follow the instructions below on how to run your ParaMonte-enabled simulations.  

+   **Running the ParaMonte simulations**  
    For complete up-to-date detailed instructions, visit: https://www.cdslab.org/paramonte/notes/run/matlab/
    +   Open the MATLAB software. On **Linux** and **macOS**, call the matlab executable from a Bash command line.  
    +   Suppose your mathematical objective function is a multivariate Normal distribution as implemented in this 
        [logfunc.m](https://raw.githubusercontent.com/cdslaborg/paramonte/master/example/mvn/MATLAB/logfunc.m) file.  
    +   For **serial** simulations, download this example generic serial 
        [main.m](https://raw.githubusercontent.com/cdslaborg/paramonte/master/example/main.m) 
        MATLAB main file and save it in the same folder containing the `logfunc.m` file that you downloaded in the above. 
        Then, simply type the name of this MATLAB main script, `main` on the MATLAB command prompt.  
    +   For **parallel** simulations, download this example generic parallel 
        [main_mpi.m](https://raw.githubusercontent.com/cdslaborg/paramonte/master/example/main_mpi.m) 
        MATLAB main file and save it in the same folder containing the `logfunc.m` file that you downloaded in the above. 
        Then, simply invoke the MPI launcher followed by the name of the MATLAB main script on a 
        MATLAB-aware MPI-aware Windows or Bash command prompt, similar to the following,  
        +   on **Windows** (preferably, on an Intel Parallel Studio command prompt or, the Microsoft Visual Studio's command prompt or, 
            some other command prompt that recognizes both `mpiexec` and `matlab` software),  
            ```  
            mpiexec -localonly -n 3 matlab -batch 'main_mpi'
            ```  
            where the `-localonly` flag is needed only if you are using the Intel MPI runtime libraries 
            (which is the default MPI library used to build the ParaMonte libraries on Windows).  
        +   on **Linux** or **macOS** (within a Bash terminal),  
            ```  
            mpiexec -n 3 matlab -batch 'main_mpi'
            ```  
        Here, the parallel simulations are performed on 3 processes. Change the number 3 to any number of processes you wish to use, 
        but do not go beyond the maximum number of physical processes available on your system, otherwise, it will only degrade 
        the performance of your parallel simulations. For example, if you are running the parallel simulation on a personal 
        quad-cores laptop, set the number of processes to either 3 or 4 at most.  
    +   Enjoy the unification of simplicity, efficiency, and parallelism in Monte Carlo simulations!  
    +   The ParaMonte library samplers are extremely versatile with many adjustable input parameters. 
        To learn about the many advanced features of the ParaMonte samplers, visit: https://www.cdslab.org/paramonte  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
