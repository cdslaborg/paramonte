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
+ On Windows and Linux operating systems, we highly recommend downloading and installing the [Intel MPI runtime libraries](https://software.intel.com/en-us/mpi-library), which is available to the public free of charge.  
+ On macOS, we recommend [Open-MPI](https://www.open-mpi.org/) since the Intel MPI library does not support macOS. For more information, visit [https://www.cdslab.org/paramonte/](https://www.cdslab.org/paramonte/).  

Citing ParaMonte  
================  

As per the license terms, we kindly ask you to cite the ParaMonte library if you use the library or any parts of it in your research, education, or software development. Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/notes/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work) for the full citation and reference information.  

Authors and contributors  
========================  

+ [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)  
    + astrophysicist/bioinformatician by training (and a science-lover in general),  
    + Ph.D. in computational physics/bioinformatics from the University of Texas at Austin,  
    + currently a faculty member of Physics and Data Science at The University of Texas at Arlington,  
    + with teaching/research experience/background in computational and data sciences, statistics, data analysis, and modeling, stochastic processes, Monte Carlo Methods, Bayesian probability theory, high energy physics, astronomy and astrophysics, computational physics, Molecular Dynamics simulations, biomedical science and MRI data analysis, bioinformatics and evolutionary biology (viral evolution, protein dynamics, and interactions),  
    + contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  

+ [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)  
    + physicist / cosmologist by training,  
    + currently a UTA Physics member,  
    + deep philosophical thinker,  
    + contact: [Fatemeh.Bagheri@uta.edu](mailto:"Fatemeh.Bagheri@uta.edu")  

License  
=======  
  
[GNU Lesser General Public License](https://github.com/cdslaborg/paramonte/blob/master/LICENSE)  
  
Example usage instructions  
==========================  

+ For complete clear organized up-to-date instructions on the build process and the installation of the ParaMonte library, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)  

## Quick start  

+ Go to the [release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases),  
+ Decide on the parallelism paradigm that you want to use: serial / MPI (the Coarray Fortran implementation is not available as a prebuilt dynamic library),  
+ Decide on the Operating System (OS) on which you want to run the ParaMonte simulations: Windows / macOS / Linux,  
+ Learn about the naming convention used for the ParaMonte prebuilt libraries [here](https://www.cdslab.org/paramonte/notes/installation/readme/#naming-convention-used-for-paramonte-library-builds),  
+ Download the prebuilt ParaMonte library of your choice based on the decisions you have made in the above,    
+ Each prebuilt library ships with a full-fledged set of example codes and build scripts. Uncompress the prebuilt library:  
    +   On **Windows**: Simply double-click on the file and select **extract files** from the Windows Explorer menu.  
    +   On **macOS/Linux**: Open a Bash terminal and navigate to the folder containing the compressed library. Use the following command to untar the compressed file,  
        ```  
        tar xvzf libparamonte
        ```  
        where you will have to replace `libparamonte` in the command above with the full name of the compressed file.  

### Building and running ParaMonte simulations on Windows  

**Environment**  

+ **Install Microsoft Visual Studio (>2017)**: You will need to have a recent Microsoft Visual Studio (MSVS) installed on your system. The community edition of this software is available free of charge. When installing MSVS, make sure to install all the C++ components of the Visual Studio.  

+ **Install Intel Parallel Studio**: If you are a student/teacher/faculty/open-source-developer, you can also download and install, free of charge, the most recent **Intel Parallel Studio** on your system which, by default, includes the Intel MPI library.  

+ **Open the right command-line interface to build/run ParaMonte example**: If the ParaMonte library that you intend to use is built for 64-bit architecture, then make sure you open a 64-bit instance of the command-line interface in either of the two cases below:  
    +   If you have installed Intel Parallel Studio, open an instance of the **command-line interface** that comes with Intel Parallel Studio from the list of programs in the Windows start menu. This is simply a Windows command prompt that has all the necessary compiler variables and paths predefined in it.  
    +   Otherwise, if you do not have Intel Parallel Studio, open an instance of the **command-line interface** that comes with Microsoft Visual Studio from the list of programs in the Windows start menu. This is simply a Windows command prompt that has all the necessary compiler variables and paths predefined in it.  

+   **Build the ParaMonte example**:  
    Build the example via Intel Parallel Studio command-line interface,  
    ```
    build.bat  
    ```
    The build script will automatically detect whether a parallel application has been built. By default, the name of the output executable is `main.exe`.  

+   **Run the ParaMonte example executable**:  
    +   For serial applications simply type the name of the output executable,  
        ```
        main.exe
        ```
    +   For parallel applications call `mpiexec`,  
        ```
        mpiexec -n NUM_PROCESSES main.exe
        ```
        where `NUM_PROCESSES` represents the number of processes on which the application will run.  

### Unix  

**Environment**  

+   If you intend to run serial ParaMonte applications, **install either Intel Fortran compiler (ifort >2018) or, GNU Fortran compiler (gfortran >7.0.0)**. If you follow the full installation instructions of ParaMonte, the prerequisites will be automatically installed for you.  

+   If you intend to run MPI parallel ParaMonte applications, **install either Intel Parallel Studio (>2018) or, GNU Compiler Collection (>7.0.0) and MPICH (>3.2) on your system**. If you follow the full installation instructions of ParaMonte, the prerequisites will be automatically installed for you. Note that on macOS, only the latter option is available since the Intel MPI library does not support macOS.  

+   If you intend to run Coarray parallel ParaMonte applications, **install either Intel Parallel Studio (>2018) or, GNU Compiler Collection (>7.0.0) and OpenCoarrays (>2.8.0) on your system**. If you follow the full installation instructions of ParaMonte, the prerequisites will be automatically installed for you. Note that on macOS, only the latter option is available since the Intel MPI library does not support macOS.  

+   Open a Bash shell, change directory to the ParaMonte example's directory, then build the executable via,  
    ```
    build.sh  
    ```
    The build script will automatically detect whether a parallel application has to be built. By default, the name of the output executable is `main.exe`. The script will also generate a new Bash script named `run.sh`. To run the generated example executable, type,  
    ```
    ./run.sh
    ```
    The script will automatically detect whether the application has to be run in parallel or serial. If the application is parallel, you can also pass the number of cores on which you want to run the example via,  
    ```
    ./run.sh --nproc NUM_PROCESSOR
    ```  
    or,  
    ```
    ./run.sh -n NUM_PROCESSOR
    ```  
    where you will have to replace `NUM_PROCESSOR` with your desired number of processes.  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
