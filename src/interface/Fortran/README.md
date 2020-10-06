
  
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

  
Installation  
============  

The ParaMonte library installation/build process is fully automated for all of the supported programming languages. 
The pre-built ready-to-use libraries are also available on [the release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases). 
Each prebuilt ParaMonte library automatically ships with a full-fledged set of example codes and build scripts.  

For more information and quick-start in the programming language of your choice, visit the [ParaMonte library homepage](https://www.cdslab.org/paramonte).  

  
Dependencies  
============  

Beyond an optional MPI runtime library for parallel simulations, the ParaMonte kernel has **zero dependency** on external third-party libraries or packages.  

  
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

+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)  
    +   physicist / cosmologist by training,  
    +   currently a UTA Physics member,  
    +   deep philosophical thinker,  
    +   contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  

  
Example usage instructions  
==========================  

+   For complete clear organized up-to-date instructions on the build process and the installation 
    of the ParaMonte library, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

## Quick start  

+   Go to the [release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases),  
+   Decide on the parallelism paradigm that you want to use: serial / MPI 
    (the Coarray Fortran implementation is not available as a prebuilt dynamic library),  
+   Decide on the Operating System (OS) on which you want to run the ParaMonte simulations: Windows / macOS / Linux,  
+   Learn about the naming convention used for the ParaMonte prebuilt libraries [here](https://www.cdslab.org/paramonte/notes/installation/readme/#naming-convention-used-for-paramonte-library-builds),  
+   Download the prebuilt ParaMonte library of your choice based on the decisions you have made in the above. 
    If you are not sure which prebuilt library is suitable for your needs, use the prebuilt library recommended 
    [here for Windows](https://www.cdslab.org/paramonte/notes/installation/windows/#using-the-prebuilt-paramonte-library), or 
    [here for Linux](https://www.cdslab.org/paramonte/notes/installation/linux/#using-the-prebuilt-paramonte-library), or 
    [here for macOS](https://www.cdslab.org/paramonte/notes/installation/macos/#using-the-prebuilt-paramonte-library).  
+   Each prebuilt library ships with a full-fledged set of example codes and build scripts. Uncompress the prebuilt library:  
    +   On **Windows**: Simply double-click on the zip-file and select **extract files** from the Windows Explorer menu.  
    +   On **macOS/Linux**: Open a Bash terminal and navigate to the folder containing the compressed library. 
        Use the following command to untar the compressed file,  
        ```  
        ls libparamonte*.tar.gz | xargs -i tar xvzf {}
        ```  
        to extract all libparamonte tar files in the current directory.  

### Building and running ParaMonte simulations on Windows  

+   **Note**: Theoretically, you can use any Fortran compiler on Windows to build and link your applications against the ParaMonte library. 
    However, the ParaMonte library example build scripts, as described below, currently only recognize the Intel Fortran compilers.  

+   **Install the Microsoft Visual Studio (>2017)**: You will need to have a recent Microsoft Visual Studio (MSVS) installed on your system. 
    The community edition of this software is available free of charge. When installing MSVS, 
    make sure to install all the C++ components of the Visual Studio.  

+   **Install the Intel Parallel Studio**: If you are a student/teacher/faculty/open-source-developer, you can also download and install, free of charge, 
    the most recent **Intel Parallel Studio** on your system which, by default, includes the Intel MPI library. You can follow the instructions given 
    [on this page](https://www.cdslab.org/recipes/programming/intel-parallel-studio-installation-windows/intel-parallel-studio-installation-windows) 
    to install the Intel Parallel Studio on your system.  

+   **Open the right command-line interface to build/run the ParaMonte example**: 
    If the ParaMonte library that you intend to use is built for 64-bit architecture, 
    then make sure you open a 64-bit instance of the command-line interface in either of the two cases below:  
    +   If you have installed Intel Parallel Studio, open an instance of the **command-line interface** that comes 
        with the Intel Parallel Studio from the list of programs in the Windows start menu. This is simply a Windows 
        command prompt that has all the necessary Intel compiler variables and paths predefined in it.  
    +   Otherwise, if you do not have Intel Parallel Studio, open an instance of the **command-line interface** that comes 
        with the Microsoft Visual Studio from the list of programs in the Windows start menu. This is simply a Windows 
        command prompt that has all the necessary compiler variables and paths predefined in it.  

+   **Build and run the ParaMonte example**:  
    Build the example via the Intel Parallel Studio command-line interface,  
    ```  
    build.bat  
    ```  
    The build script will automatically detect whether a parallel simulation has been built. 
    By default, the name of the output executable is `main.exe`. 
    *Note that the build script will both build and run the executable*.  

+   **Run the ParaMonte example executable**:  
    +   For serial simulations, simply type the name of the output executable,  
        ```  
        main.exe
        ```  
    +   For parallel simulations, invoke the MPI launcher `mpiexec`,  
        ```  
        mpiexec -n NUM_PROCESSES main.exe
        ```  
        where `NUM_PROCESSES` represents the number of processes on which the simulation will run. 
        If you are using the Intel MPI library to run your ParaMonte application in parallel, we also recommend using the `-localonly` flag. 
        See [this page](https://www.cdslab.org/paramonte/notes/run/#running-the-manually-generated-executable-on-multiple-processors-on-windows) 
        for usage and utilities of this Intel MPI launcher flag.  

### Building and running ParaMonte simulations on macOS / Linux  

+   **Note**: Theoretically, you can use any Fortran compiler on macOS/Linux to build and link your applications against the ParaMonte library. 
    However, the ParaMonte library example build scripts, as described below, currently only recognize the Intel and GNU Fortran compilers.  

+   If you intend to run **serial** ParaMonte simulations, install either,  
    +   **the Intel Fortran compiler (ifort >2018)**, or,  
    +   **the GNU Fortran compiler (gfortran >7.0.0)**,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  

+   If you intend to run **MPI parallel** ParaMonte simulations, install either,  
    +   **the Intel Parallel Studio (>2018)** on Linux, or,  
    +   **the GNU Compiler Collection (>7.0.0) and the MPICH (>3.2) library** on Linux,  
    +   **the GNU Compiler Collection (>7.0.0) and the OpenMPI (>4.0) library** on macOS,  
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  
    Note that on **macOS**, only the latter option (the GNU compilers) 
    is available since the Intel MPI library does not support the macOS platform.  

+   If you intend to run Coarray parallel ParaMonte simulations, install either,  
    +   **the Intel Parallel Studio (>2018)**, or,  
    +   **the GNU Compiler Collection (>7.0.0) and OpenCoarrays (>2.8.0)**, 
    on your system. If you follow the full installation instructions of the ParaMonte library, 
    these compiler prerequisites will be automatically installed for you.  
    Note that on **macOS**, only the latter option (the GNU compilers) 
    is available since the Intel MPI library does not support the macOS platform.  

+   Open a Bash terminal, change directory to the ParaMonte prebuilt library's directory, then build the executable via,  
    ```  
    build.sh  
    ```  
    The build script will automatically detect whether a parallel simulation has to be built. 
    By default, the name of the output executable is `main.exe`. The script will also generate 
    a new Bash script named `run.sh`. To run the generated example executable, type,  
    ```  
    ./run.sh
    ```  
    The script will automatically detect whether the simulation has to be run in parallel or serial. 
    If the simulation is parallel, you can also pass the number of cores on which you want to run the example via,  
    ```  
    ./run.sh --nproc NUM_PROCESSOR
    ```  
    or,  
    ```  
    ./run.sh -n NUM_PROCESSOR
    ```  
    where you will have to replace `NUM_PROCESSOR` with your desired number of processes.  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
