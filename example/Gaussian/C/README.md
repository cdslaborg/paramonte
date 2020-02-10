> **ParaMonte: plain powerful parallel Monte Carlo library.**  
> 
> Copyright (C) 2012-present, The Computational Data Science Lab  
> 
> This file is part of the ParaMonte library.   
> 
> ParaMonte is free software: you can redistribute it and/or modify  
> it under the terms of the GNU Lesser General Public License as published by  
> the Free Software Foundation, version 3 of the License.  
> 
> ParaMonte is distributed in the hope that it will be useful,  
> but WITHOUT ANY WARRANTY; without even the implied warranty of  
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  
> GNU Lesser General Public License for more details.  
> 
> You should have received a copy of the GNU Lesser General Public License  
> along with ParaMonte.  If not, see [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).  
  

**NOTE:** For complete clear instructions on the installation and building of ParaMonte, please visit [cdslab.org/pm](https://www.cdslab.org/pm). All of the steps described below will be automatically done for you if you build the entire library by following the instructions provided in the aforementioned web-link.  

## Build  

Building this example requires the following components in the same folder and the following software installed on your system:  

- the example source files,  
- the example input file `paramonte.in`,  
- on Windows,  
    - the `build.bat` script, available in the grandparent directory of this folder,  
    - either Intel C/C++ (`icl`) or Microsoft Visual C++ (`cl`) compiler,  
    - optionally, if running in parallel, the installation of Intel MPI library is required, which is available free of charge,  
- on Linux,  
    - either Intel C/C++ (`icc/icpc`) or GNU C/C++ (`gcc/g++`) compiler,  
    - the `build.sh` script, available in the grandparent directory of this folder,  
    - optionally, if running in parallel, the installation of Intel MPI library is required, which is available free of charge,  
- on Mac,  
    - GNU C/C++ (`gcc/g++`) compiler,  
    - the `build.sh` script, available in the grandparent directory of this folder,  
    - optionally, if running in parallel, the installation of OpenMPI library is required, which is available free of charge,  

## Usage  

### Windows  

**Environment**  

- **Install Microsoft Visual Studio (>2017)**: You will need to have a recent Microsoft Visual Studio (MSVS) installed on your system. The community edition of this software is available free of charge. When installing MSVS, make sure to install all the C++ components of the Visual Studio.  

- **Optionally, install Intel Parallel Studio**: If you are a student/teacher/faculty/open-source-developer, you can also download and install, free of charge, the most recent **Intel Parallel Studio** on your system which, by default, includes the Intel MPI library.  

- **Open the right command-line interface to build/run ParaMonte example**: If the ParaMonte library that you intend to use is built for 64-bit architecture, then make sure you open a 64-bit instance of the command-line interface in either of the two cases below:  
    - If you have installed Intel Parallel Studio, open an instance of the **command-line interface** that comes with Intel Parallel Studio from the list of programs in the Windows start menu. This is simply a Windows command prompt that has all the necessary compiler variables and paths predefined in it.  
    - Otherwise, if you do not have Intel Parallel Studio, open an instance of the **command-line interface** that comes with Microsoft Visual Studio from the list of programs in the Windows start menu. This is simply a Windows command prompt that has all the necessary compiler variables and paths predefined in it.  

- **Build the ParaMonte example**:  
    - To build the example via Intel Parallel Studio,  
        ```
        build.bat  
        ```
    - To build the example via Microsoft Visual Studio,  
        ```
        build.bat msvc  
        ```
        where the passed argument `msvc` implies the use of Microsoft Visual C++ compiler for building the application. The build script will automatically detect whether a parallel application has been built. By default, the name of the output executable is `runExample.exe`.  

- **Run the ParaMonte example executable**:  
    - For serial applications simply type the name of the output executable,  
        ```
        runExample.exe
        ```
    - For parallel applications call `mpiexec`,  
        ```
        mpiexec -np NUM_PROCESSES runExample.exe
        ```
        where `NUM_PROCESSES` represents the number of processes on which the application will run.  

### Unix  

**Environment**  

-   If you intend to run serial ParaMonte applications, **install either Intel C/C++ compiler (icc/icpc >2018) or, GNU Fortran compiler (gcc/g++ >7.0.0)**. If you follow the full installation instructions of ParaMonte, the prerequisites will be automatically installed for you.  

-   If you intend to run MPI parallel ParaMonte applications, **install either Intel Parallel Studio (>2018) or, GNU Compiler Collection (>7.0.0) and MPICH (>3.2) on your system**. If you follow the full installation instructions of ParaMonte, the prerequisites will be automatically installed for you. Note that on macOS, only the latter option is available since the Intel MPI library does not support macOS.  

-   Open a Bash shell, change directory to the ParaMonte example's directory, then build the executable via,  
    ```
    build.sh  
    ```
    The build script will automatically detect whether a parallel application has to be built. By default, the name of the output executable is `runExample.exe`. The script will also generate a new Bash script named `run.sh`. To run the generated example executable, type,  
    ```
    ./run.sh
    ```
    The script will automatically detect whether the application has to be run in parallel or serial. If the application is parallel, you can also pass the number of cores on which you want to run the example via,  
    ```
    ./run.sh --num_images NUM_PROCESSOR
    ```  
    or,  
    ```
    ./run.sh -n NUM_PROCESSOR
    ```  
    where you will have to replace `NUM_PROCESSOR` with your desired number of processes.  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/pm) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
