#!/usr/bin/env python
####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################

#### ATTN: This code must be executed from the base directory where it exists. 
#### Navigate the directory containing this file, then: python getReadme.py

from getShield import getShield
banner = getShield()
#with open("shields.html", "r") as file: banner = file.read()

banner = """

<div align="center">
<a href="https://www.cdslab.org/paramonte" target="_blank">
    <img src="https://github.com/cdslaborg/paramonte/blob/main/img/paramonte.png" alt="ParaMonte: Parallel Monte Carlo and Machine Learning Library" style="width:500px;height:500px;" />
</a>
<br>
""" + banner + """
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
"""

####################################################################################################################################

class Struct: pass

sectionList =   [ "quickstart"
                , "installation"
                , "dependencies"
                , "parallelism"
                , "examples"
                , "citation"
                , "license"
                , "authors"
                , "contribution"
                ]

readme = dict()
for section in sectionList: readme[section] = dict()

####################################################################################################################################

readme["installation"]["title"] = """

Installation  
============  
"""

readme["installation"]["main"] = readme["installation"]["title"] + """

The prebuilt, ready-to-use libraries are available on 
[the release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases).  
Each prebuilt ParaMonte library automatically ships with a full-fledged set of example codes and build scripts.  

Alternatively, you can build the library from the source 
in the project's [GitHub repository](https://github.com/cdslaborg/paramonte).  
The ParaMonte library installation/build process is fully automated for all supported programming languages.  
Currently, the following compiler suites are supported **for builds from source**,  

| Compiler Suite                    | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| GNU Compiler Collection >= 10.3   | ✅              | ✅            | ✅             | ✅            |  
| Intel OneAPI >= 2021.8.0          | ✅              | ✅            | ✅             | ❌            |  

with the following essential build systems,

| Dependency                        | Windows (amd64) | Linux (amd64) | macOS (amd64) | macOS (arm64) |  
|----------------------------------:|:---------------:|:-------------:|:-------------:|:-------------:|  
| CMake >= 3.21                     | ✅              | ✅            | ✅             | ✅            |  

+   Read the [`install.md` installation instruction notes](install.md) 
    to learn about the library's optional build configurations.  
+   For a quickstart in the programming language of your choice, 
    visit the [ParaMonte library homepage](https://www.cdslab.org/paramonte).  
"""

readme["installation"]["c"] = readme["installation"]["main"]

readme["installation"]["cpp"] = readme["installation"]["main"]

readme["installation"]["fortran"] = readme["installation"]["main"]

readme["installation"]["matlab"] = readme["installation"]["title"] + """

We strongly advise you to follow the guidelines in 
[this documentation page](https://www.cdslab.org/paramonte/generic/1/installation/matlab/) 
to either download the pre-built ParaMonte MATLAB libraries for your platform or to build the library from the source code.  
"""

readme["installation"]["python"] = readme["installation"]["title"] + """
The latest release of ParaMonte can be installed from PyPI using `pip`:

    pip3 install --user --upgrade paramonte

or,

    pip install --user --upgrade paramonte

Alternatively, you can follow the guidelines in 
[this documentation page](https://www.cdslab.org/paramonte/generic/1/installation/python/) 
to either download the pre-built ParaMonte Python libraries for your platform or to build the library from the source code.  
"""

####################################################################################################################################

readme["dependencies"]["title"] = """

Dependencies
============

Aside from an optional MPI runtime library for MPI-parallel simulations, 
the ParaMonte library core has **zero dependency** on external third-party libraries or packages.  
"""

readme["dependencies"]["main"] = readme["dependencies"]["title"] + """
"""

readme["dependencies"]["c"] = readme["dependencies"]["main"]

readme["dependencies"]["cpp"] = readme["dependencies"]["main"]

readme["dependencies"]["fortran"] = readme["dependencies"]["main"]

readme["dependencies"]["matlab"] = readme["dependencies"]["title"] + """
"""

readme["dependencies"]["python"] = readme["dependencies"]["title"] + """
The Python interface of ParaMonte depends on very few third-party libraries. 
These include `numpy`, `scipy`, `pandas`, and `matplotlib`. 
The last two (plotting) libraries are only used for the post-processing of simulation 
results and are therefore not needed if you do not plan to use the post-processing 
features of the ParaMonte library. If you have a recent version of the 
[Anaconda Python ](https://www.anaconda.com/download) distribution 
installed on your system, then all the dependencies already 
exist and are automatically installed.  
"""

####################################################################################################################################

readme["parallelism"]["title"] = """

Parallelism
===========
"""

readme["parallelism"]["main"] = readme["parallelism"]["title"] + """
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
"""

readme["parallelism"]["c"] = readme["parallelism"]["main"]

readme["parallelism"]["cpp"] = readme["parallelism"]["main"]

readme["parallelism"]["fortran"] = readme["parallelism"]["main"]

readme["parallelism"]["matlab"] = readme["parallelism"]["main"] + """

+   The serial version of the ParaMonte MATLAB library has **NO external library dependencies**.  
+   The MPI-parallel version requires a ParaMonte-compatible MPI runtime library.  
+   The new thread-based parallelism features of the ParaMonte MATLAB library require the MATLAB Parallel Toolbox.  
+   The MATLAB Parallel Toolbox is usually automatically shipped with MATLAB software.  
"""

readme["parallelism"]["python"] = readme["parallelism"]["main"]

####################################################################################################################################

readme["quickstart"]["title"] = """

Quick Start
===========
"""

readme["quickstart"]["main"] = readme["quickstart"]["title"] + """
+   Follow the quick start instructions in this [QUICKSTART.md](https://github.com/cdslaborg/paramonte/blob/main/QUICKSTART.md) file. 
+   For programming languages other than C, C++, Fortran, see the quick start section of the README.md file in the corresponding 
    language's source subfolder in the [src folder](https://github.com/cdslaborg/paramonte/src) in the root directory of the ParaMonte repository.
"""

readme["quickstart"]["c"] = readme["quickstart"]["main"]

readme["quickstart"]["cpp"] = readme["quickstart"]["main"]

readme["quickstart"]["fortran"] = readme["quickstart"]["main"]

#> The MATLAB branch of the latest release of ParaMonte is a work in progress.  
#> While most library functionalities are already implemented, the library API remains unstable and subject to changes.  
#> We encourage you to build the library from source on your system and test or use it and let us know any comments you may have.  
#> Remember that the previous major release of ParaMonte MATLAB remains fully functional, albeit with entirely different API.  
#> The following guidelines apply to the previous major release of [ParaMonte MATLAB](https://www.cdslab.org/paramonte/generic/1/installation/matlab/).  
#The corresponding example source files (the `*.mlx` files) can be downloaded from the [paramontex GitHub repository](https://github.com/cdslaborg/paramontex/tree/main/MATLAB/mlx),
#a repository dedicated to the ParaMonte library examples.
readme["quickstart"]["matlab"] = readme["quickstart"]["title"] + """

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
"""

readme["quickstart"]["python"] = readme["quickstart"]["title"] + """
For a quick start with some Jupyter Notebook examples, visit [this ParaMonte documentation page](https://www.cdslab.org/paramonte/generic/latest/examples/python/jupyter/).
The corresponding example source files (the `*.ipynb` files) can be downloaded from the [paramontex GitHub repository](https://github.com/cdslaborg/paramontex/tree/main/Python/Jupyter),
which is a repository dedicated to the ParaMonte library examples.

The following example code samples a 4-dimensional MultiVariate Normal (MNV) distribution via the ParaDRAM sampler in serial mode,

```python
import numpy as np
import paramonte as pm
def getLogFunc(point): return -0.5 * np.dot(point, point)
sampler = pm.sampling.Paradram()
sampler.run ( getLogFunc = getLogFunc   # the objective function
            , ndim = 4 # assume 4-dimensional objective function
            )
```

To learn about the ParaMonte Python library's post-processing and visualization tools, visit [this documentation page](https://www.cdslab.org/paramonte/generic/latest/examples/python/jupyter/).
"""

####################################################################################################################################

readme["citation"]["title"] = """

Citing ParaMonte
================
"""

readme["citation"]["main"] = readme["citation"]["title"] + """
The ParaMonte library is an honor-ware and its currency is acknowledgment and citations. 

If you use ParaMonte or any ideas from this software, please acknowledge it by citing the ParaMonte library's 
main publications as listed in [ACKNOWLEDGMENT.md](https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md).  

Visit [the ParaMonte library homepage](https://www.cdslab.org/paramonte/generic/latest/overview/preface/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work)
to access the PDF version of these files free of charge.  
"""

readme["citation"]["c"] = readme["citation"]["main"]

readme["citation"]["cpp"] = readme["citation"]["main"]

readme["citation"]["fortran"] = readme["citation"]["main"]

readme["citation"]["matlab"] = readme["citation"]["main"]

readme["citation"]["python"] = readme["citation"]["main"]

####################################################################################################################################

readme["license"]["title"] = """

License
=======
"""

readme["license"]["main"] = readme["license"]["title"] + """
The majority of the ParaMonte library routines are distributed under the permissive 
[MIT License](LICENSE.md) plus acknowledgments that are detailed in [LICENSE.md](./LICENSE.md).  

However, there are occasionally modernized and extended routines from external projects, particularly within the ParaMonte Fortran library. 
These routines may have a license that does not fully overlap with the default ParaMonte library license [LICENSE.md](./LICENSE.md). 
In such cases, the documentation and source code of the specific routine explicitly list the original license information. 

This is free software, so help us keep it freely available to the public by acknowledging the usage and contributing to it. 
If you have questions or concerns about the license, please get in touch with the project's lead (shahmoradi@utexas.edu). 
"""

readme["license"]["c"] = readme["license"]["main"]

readme["license"]["cpp"] = readme["license"]["main"]

readme["license"]["fortran"] = readme["license"]["main"]

readme["license"]["matlab"] = readme["license"]["main"]

readme["license"]["python"] = readme["license"]["main"]

#################################################################################################################################################################

readme["authors"]["title"] = """

Authors and Contributors
========================
"""

amir = """
+   [Amir Shahmoradi](https://www.cdslab.org/people/#amir-shahmoradi)
    +   Astrophysicist/bioinformatician by training,
    +   Ph.D. in computational physics/bioinformatics, The University of Texas at Austin,
    +   Currently a faculty member of the Data Science and Physics programs at UT Arlington,
    +   Contact: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")
"""

fatima = """
+   [Fatemeh Bagheri](https://www.linkedin.com/in/fbagheri)
    +   PhD in Astronomy,
    +   PhD in Space Physics,
    +   Deep philosophical thinker,
    +   Currently at NASA Goddard Space Flight Center,
    +   Contact: [bagheri.fateme@gmail.com](mailto:"bagheri.fateme@gmail.com")
"""

josh = """
+   [Joshua Osborne](https://www.cdslab.org/people/#joshua-alexander-osborne)
    +   Physicist / computational data scientist by training,
    +   Currently a Physics PhD candidate at UT Arlington,
    +   Contact: [joshuaalexanderosborne@gmail.com](mailto:"joshuaalexanderosborne@gmail.com")
"""

readme["authors"]["main"] = readme["authors"]["title"] + amir + fatima + josh

readme["authors"]["c"] = readme["authors"]["title"] + amir

readme["authors"]["cpp"] = readme["authors"]["title"] + amir + fatima

readme["authors"]["fortran"] = readme["authors"]["title"] + amir + fatima

readme["authors"]["matlab"] = readme["authors"]["title"] + amir + josh + fatima

readme["authors"]["python"] = readme["authors"]["title"] + fatima + amir + josh

####################################################################################################################################

readme["examples"]["title"] = """

Example Usage Instructions
==========================
"""

readme["examples"]["main"] = readme["examples"]["title"] + """
+   The best way to get started with example usage is to build the language-specific library examples.  
    The library example build and run scripts are automatically generated for each ParaMonte build and 
    are automatically stored in the output installation directory of the build.

+   For complete, organized, up-to-date instructions, visit: [cdslab.org/pm](https://www.cdslab.org/paramonte)

+   For a quick look into *language-specific* README.md instructions, visit:
    +   **C**: [https://github.com/cdslaborg/paramonte/tree/main/src/c](https://github.com/cdslaborg/paramonte/tree/main/src/c)
    +   **C++**: [https://github.com/cdslaborg/paramonte/tree/main/src/cpp](https://github.com/cdslaborg/paramonte/tree/main/src/cpp)
    +   **Fortran**: [https://github.com/cdslaborg/paramonte/tree/main/src/fortran](https://github.com/cdslaborg/paramonte/tree/main/src/fortran)
    +   **MATLAB**: [https://github.com/cdslaborg/paramonte/tree/main/src/matlab](https://github.com/cdslaborg/paramonte/tree/main/src/matlab)

+   The ParaMonte library for Python is currently undergoing internal testing before their impending release. 
    For now, the previous versions (1) of the ParaMonte Python libraries can be downloaded and used from the previous 
    library major release on [GitHub Release Page](https://github.com/cdslaborg/paramonte/releases/tag/v1.5.1).
"""

### Quick start
#
#+   Go to the [release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases).
#+   Decide on the parallelism paradigm you want to use: serial / MPI / OpenMP.
#+   Decide on the operating system (OS) you want to run the ParaMonte simulations: Windows, macOS, or Linux.
#+   Learn about the naming convention used for the ParaMonte prebuilt libraries
#    [here](https://www.cdslab.org/paramonte/generic/latest/installation/readme/#naming-convention-used-for-paramonte-library-builds).
#+   Download the prebuilt ParaMonte library of your choice based on the decisions you have made above.
#    If you are unsure which prebuilt library is suitable for your needs, use the prebuilt library recommended
#    [here for Windows](https://www.cdslab.org/paramonte/generic/latest/installation/windows/#using-the-prebuilt-paramonte-library),
#    or [here for Linux](https://www.cdslab.org/paramonte/generic/latest/installation/linux/#using-the-prebuilt-paramonte-library),
#    or [here for macOS](https://www.cdslab.org/paramonte/generic/latest/installation/macos/#using-the-prebuilt-paramonte-library).
#+   Each prebuilt library ships with a full-fledged set of example codes and build scripts. Uncompress the prebuilt library:
#    +   On **Windows**: Simply double-click on the zip-file and select **extract files** from the Windows Explorer menu.
#    +   On **macOS/Linux**: Open a Bash terminal and navigate to the compressed library folder.
#        Use the following command to untar the compressed file,
#        ```
#        ls libparamonte*.tar.gz* | xargs -i tar xvzf {}
#        ```
#        to extract all `libparamonte` tar files in the current directory.
#
#### Building and running ParaMonte simulations  
#
#+   **Note**: Theoretically, you can use any C/C++ compiler on Windows to build and link your applications against the ParaMonte library. 
#    However, remember that the ParaMonte library example build scripts have only been tested with the GNU, Microsoft, and Intel C/C++ compilers. 
#
#+   **Note**: Theoretically, you can use any C/C++ compiler on macOS or Linux to build and link your applications against the ParaMonte library. 
#    However, as described below, the ParaMonte library example build scripts only recognize the Intel and GNU C/C++ compilers. 
#
#+   For the generic library build instructions, follow the guidelines provided in [install.md](../../install.md) in the root directory of the repository. 
readme["examples"]["c"] = readme["examples"]["title"] + """
+   For complete clear organized, up-to-date instructions on the build process and the 
    installation of the ParaMonte library, visit the [ParaMonte C examples guidelines](https://www.cdslab.org/paramonte/generic/latest/usage/examples/c/).  
"""

readme["examples"]["cpp"] = readme["examples"]["title"] + """
+   For complete clear organized, up-to-date instructions on the build process and the 
    installation of the ParaMonte library, visit the [ParaMonte C++ examples guidelines](https://www.cdslab.org/paramonte/generic/latest/usage/examples/cpp/).  
"""

### Quick start
#
#+   Go to the [release page of the ParaMonte library on GitHub](https://github.com/cdslaborg/paramonte/releases),
#+   Decide on the parallelism paradigm that you want to use: serial / MPI / OpenMP
#    (the Coarray Fortran implementation is not available as a prebuilt dynamic library),
#+   Decide on the Operating System (OS) on which you want to run the ParaMonte simulations: Windows / macOS / Linux,
#+   Learn about the naming convention used for the ParaMonte prebuilt libraries 
#    [here](https://www.cdslab.org/paramonte/generic/latest/installation/readme/#naming-convention-used-for-paramonte-library-builds),
#+   Download the prebuilt ParaMonte library of your choice based on the decisions you have made above.
#    If you are unsure which prebuilt library is suitable for your needs, use the prebuilt library recommended
#    [here for Windows](https://www.cdslab.org/paramonte/generic/latest/installation/windows/#using-the-prebuilt-paramonte-library), or
#    [here for Linux](https://www.cdslab.org/paramonte/generic/latest/installation/linux/#using-the-prebuilt-paramonte-library), or
#    [here for macOS](https://www.cdslab.org/paramonte/generic/latest/installation/macos/#using-the-prebuilt-paramonte-library).
#+   Each prebuilt library ships with a full-fledged set of example codes and build scripts. Uncompress the prebuilt library:
#    +   On **Windows**: Simply double-click on the zip-file and select **extract files** from the Windows Explorer menu.
#    +   On **macOS/Linux**: Open a Bash terminal and navigate to the compressed library folder.
#        Use the following command to untar the compressed file,
#        ```
#        ls libparamonte*.tar.gz | xargs -i tar xvzf {}
#        ```
#        to extract all `libparamonte` tar files in the current directory.
#
#### Building and running ParaMonte simulations  
#
#+   **Note**: Theoretically, you can use any Fortran compiler on Windows, Linux, 
#    and macOS to build and link your applications against the ParaMonte library.
#    However, remember that the ParaMonte library example build scripts 
#    have only been tested with the GNU and Intel Fortran compilers.  
#
#+   For the generic library build instructions, follow the guidelines provided 
#    in [install.md](../../install.md) in the root directory of the repository.
readme["examples"]["fortran"] = readme["examples"]["title"] + """
+   For complete clear organized, up-to-date instructions on the build process and the 
    installation of the ParaMonte library, visit the [ParaMonte Fortran examples guidelines](https://www.cdslab.org/paramonte/generic/latest/usage/examples/fortran/).  
"""

#
#+   **Install a MATLAB >2017b distribution**, preferably, the latest MATLAB.
#    The ParaMonte MATLAB library has been tested only with MATLAB version 2018b and newer.
#
#+   **Optionally install a compatible MPI library** (or let the ParaMonte library take care of the MPI installation
#    when you call the library for the first time). For parallel simulations (via MPI), you will need an MPI library
#    already installed on your system. If you choose to install the library by yourself, we recommend the Intel MPI
#    library, which is available for free from the Intel website or from [the ParaMonte GitHub release page](https://github.com/cdslaborg/paramonte/releases/tag/v1.5.1).
#    On macOS, the OpenMPI (or MPICH) MPI library can be used instead of the Intel MPI library, which currently does not support macOS.
#
#+   **Calling the ParaMonte library for the first time**
#    +   Depending on your platform,
#        +   **Windows**
#            Nothing special needs to be done. You are all set! To call the ParaMonte library for the first time, follow the instructions below.
#        +   **Linux/macOS**
#            Open a `Bash` or `zsh` terminal and open MATLAB from the command line by calling its name,
#            ```bash
#            matlab
#            ```
#            If `matlab` is not recognized on your command line as an application, seek help from
#            [this ParaMonte documentation page](https://www.cdslab.org/paramonte/generic/latest/troubleshooting/bash-matlab-command-not-found/).
#    +   Once the MATLAB interactive environment opens, navigate to the root folder of the ParaMonte library (where the LICENSE file exists)
#        and call the ParaMonte library for the first time via the following commands (simply type the commands on the MATLAB command prompt),
#        ```matlab
#        addpath(genpath("./"),"-begin"); % add the ParaMonte library directories to MATLAB's list of search paths.
#        pm = paramonte(); % instantiate an object of class paramonte.
#        pm.verify(); % verify the integrity of the ParaMonte library on your system.
#        ```
#        If needed, follow any extra instructions provided by the library on your MATLAB command prompt.
#        **If you do not intend to run simulations in parallel, you can say NO (`n`) to any MPI library installation requests via ParaMonte**.
#        If you do not intend to run simulations in parallel, answer YES (`y`) to any permission requests by
#        the ParaMonte library to install the MPI libraries on your system.
#+   **Running the ParaMonte simulations**
#    For complete, up-to-date, detailed instructions, visit: https://www.cdslab.org/paramonte/generic/latest/run/matlab/
#    +   Open the MATLAB software. On **Linux** and **macOS**, call the matlab executable from a Bash command line.
#    +   The ParaMonte library typically ships with example scripts. If you see a file named `main.m` at the root directory of
#        your ParaMonte library, then call this MATLAB script to run the example simulation provided by the library.
#        If no such file exists, or if you intend to do simulations in parallel, then follow the rest of the instructions below.
#    +   Suppose your mathematical objective function is a multivariate Normal distribution as implemented in this
#        [logfunc.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/mvn/MATLAB/logfunc.m) file.
#    +   For **serial** simulations, download this example generic serial
#        [main.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main.m)
#        MATLAB main file and save it in the same folder containing the `logfunc.m` file you downloaded above.
#        Then, type the name of this MATLAB main script, `main` on the MATLAB command prompt.
#    +   For **parallel** simulations, download this example generic parallel
#        [main_mpi.m](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main_mpi.m) MATLAB main file and save it in the same folder containing the `logfunc.m` file that you downloaded above.
#        Then, invoke the MPI launcher followed by the name of the MATLAB main script on a
#        MATLAB-aware MPI-aware Windows or Bash command prompt, similar to the following,
#        +   on **Windows** (preferably, on an Intel Parallel Studio command prompt or, the Microsoft Visual Studio's command prompt or,
#            some other command prompt that recognizes both `matlab` and the Intel's `mpiexec` software),
#            ```
#            mpiexec -localonly -n 3 matlab -batch main_mpi
#            ```
#            where the `-localonly` flag is needed only if you are using the Intel MPI runtime libraries
#            (which is the default MPI library used to build the ParaMonte libraries on Windows).
#        +   on **Linux** or **macOS** (within a Bash terminal),
#            ```
#            mpiexec -n 3 matlab -batch main_mpi
#            ```
#        Here, the parallel simulations are performed on three processes. Change the number 3 to any number of processes you wish to use,
#        but do not exceed the maximum number of physical processes available on your system. Otherwise, it will only degrade
#        the performance of your parallel simulations. For example, if you run the parallel simulation on a personal
#        quad-core laptop, set the number of processes to 3 or 4 at most.
#    +   Enjoy the unification of simplicity, efficiency, and parallelism in Monte Carlo simulations!
#    +   The ParaMonte library samplers are extremely versatile, with many adjustable input parameters.
#        To learn about the many advanced features of the ParaMonte samplers, visit: https://www.cdslab.org/paramonte
readme["examples"]["matlab"] = readme["examples"]["title"] + """
+   For complete clear organized, up-to-date instructions on the build process and the 
    installation of the ParaMonte library, visit the [ParaMonte MATLAB examples guidelines](https://www.cdslab.org/paramonte/generic/latest/usage/examples/matlab/).  
"""

#+   **Install a Python 3 distribution**, preferably, the Anaconda distribution of Python.
#    The Anaconda distribution of Python automatically ships with all of the
#    ParaMonte Python package dependencies when installed on your system.
#
#+   **Optionally install a compatible MPI library** (or let the ParaMonte library take care of the installation
#    when you import the package into your Python session for the first time). For parallel simulations (via MPI),
#    you will need an MPI library already installed on your system. If you choose to build the library yourself,
#    we recommend the Intel MPI library, free from the Intel website. The OpenMPI
#    library can be used on macOS instead of the Intel MPI library, which currently does not support macOS.
#
#+   **Running the ParaMonte simulations**
#    +   Open an Anaconda command-line interface or `jupyter` notebook.
#    +   Suppose your mathematical objective function is a multivariate Normal distribution as implemented in this
#        [logfunc.py](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/mvn/Python/logfunc.py) file.
#    +   For **serial** simulations, download this example generic serial
#        [main.py](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main.py) Python main file and save it in the same folder containing the `logfunc.py` file that you downloaded above.
#        Then, type the name of the Python main script, `python main.py` on the Bash terminal or the Anaconda command line.
#    +   For **parallel** simulations, download this example generic parallel
#        [main_mpi.py](https://raw.githubusercontent.com/cdslaborg/paramonte/main/example/main_mpi.py) Python main file and save it in the same folder containing the `logfunc.py` file that you downloaded in the above. Then, invoke the MPI launcher followed by the name of the Python main script on the Bash terminal, similar to the following,
#        +   on Windows (within the Anaconda command line or a terminal that recognizes both `mpiexec` and `python` software),
#            ```
#            mpiexec -localonly -n 3 python main_mpi.py
#            ```
#            where the `-localonly` flag is needed only if you are using the Intel MPI runtime libraries (which is the default MPI library used to build the ParaMonte libraries on Windows).
#        +   on macOS or Linux (within a Bash terminal),
#            ```
#            mpiexec -n 3 python main_mpi.py
#            ```
#        Here, the parallel simulations are performed on three processes. Change the number 3 to any number of processes you wish to use,
#        but do not exceed the maximum number of physical processes available on your system. Otherwise, it will only degrade
#        the performance of your parallel simulations. For example, if you run the parallel simulation on a personal
#        quad-core laptop, set the number of processes to 3 or 4 at most.
#    +   Enjoy combining simplicity, efficiency, and parallelism in Monte Carlo simulations!
#    +   The ParaMonte library samplers are incredibly versatile, with many adjustable input parameters.
#        To learn about the many advanced features of the ParaMonte routines, visit: https://www.cdslab.org/paramonte
readme["examples"]["python"] = readme["examples"]["title"] + """
+   For complete clear organized, up-to-date instructions on the build process and the 
    installation of the ParaMonte library, visit the [ParaMonte Python examples guidelines](https://www.cdslab.org/paramonte/generic/latest/usage/examples/python/).  
"""

####################################################################################################################################

readme["contribution"]["title"] = """

Contributing to ParaMonte
=========================
"""

readme["contribution"]["main"] = readme["contribution"]["title"] + """
Please read the generic contribution guidelines listed in [CONTRIBUTING.md](./CONTRIBUTING.md).
"""

readme["contribution"]["c"] = readme["contribution"]["main"]

readme["contribution"]["cpp"] = readme["contribution"]["main"]

readme["contribution"]["fortran"] = readme["contribution"]["main"]

readme["contribution"]["matlab"] = readme["contribution"]["main"]

readme["contribution"]["python"] = readme["contribution"]["main"]

####################################################################################################################################

finalNote = """

**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/paramonte) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")
"""

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

import os

abspath = os.path.abspath(__file__)
cwd = os.path.dirname(abspath)
os.chdir(cwd)
print(cwd)

fileDirDict =   { "main"    : ".."
                , "c"       : "../src/c"
                , "cpp"     : "../src/cpp"
                , "fortran" : "../src/fortran"
                , "matlab"  : "../src/matlab"
                , "python"  : "../src/python"
                }

for lang, fileDir in fileDirDict.items():

    print("Generating README.md file for the {} page...".format(lang))

    if os.path.isdir(fileDir):
        filePath = os.path.join(fileDir, "README.md")
        with open(filePath, "w") as file:
            thisReadme = ""
            for section in sectionList:
                print("    " + section)
                thisReadme += readme[section][lang]
            file.write("{}".format(banner + thisReadme + finalNote))

####################################################################################################################################
