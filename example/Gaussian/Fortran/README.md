> **ParaMonte: plain powerful parallel Monte Carlo library.**  
> 
> Copyright (C) 2012-present, The Computational Data Science Lab  
> 
> This file is part of ParaMonte library.   
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
  

**NOTE:** For complete clear instructions on installation and building of ParaMonte, please visit [cdslab.org/pm](https://www.cdslab.org/pm). All of the steps described below will be automatically done for you if you build the entire library by following the instructions provided in the aforementioned web-link.  

## Build  

Building this example requires the following components in the same folder and the following software installed on your system:  

- the example source files,  
- the example input file `paramonte.in`,  
- on Windows,  
    - the `build.bat` script, available in the grandparent directory of this folder,  
    - currently the only compiler that the build script recognizes is Intel Fortran (`ifort`) compiler,  
    - optionally, if running in parallel, the installation of Intel MPI library is required, which is available free of charge,  
- on Linux,  
    - either Intel Fortran (`ifort`) or GNU Fortran (`gfortran`) compiler,  
    - the `build.sh` script, available in the grandparent directory of this folder,  
    - optionally, if running in parallel, the installation of Intel MPI library is required, which is available free of charge,  
- on Mac,  
    - GNU Fortran (`gfortran`) compiler,  
    - the `build.sh` script, available in the grandparent directory of this folder,  
    - optionally, if running in parallel, the installation of OpenMPI library is required, which is available free of charge,  

## Usage  

### Windows  

**Environment**  

- **Install Microsoft Visual Studio (>2017)**: You will need to have a recent Microsoft Visual Studio (MSVS) installed on your system. The community edition of this software is available free of charge. When installing MSVS, make sure to install all the C++ components of the Visual Studio.  

- **Install Intel Parallel Studio**: If you are a student/teacher/faculty/open-source-developer, you can also download and install, free of charge, the most recent **Intel Parallel Studio** on your system which, by default, includes the Intel MPI library.  

- **Open the right command-line interface to build/run ParaMonte example**: If the ParaMonte library that you intend to use is built for 64-bit architecture, then make sure you open a 64-bit instance of the command-line interface in either of the two cases below:  
    - If you have installed Intel Parallel Studio, open an instance of the **command-line interface** that comes with Intel Parallel Studio from the list of programs in Windows start menu. This is simply a Windows command prompt that has all the necessary compiler variables and paths predefined in it.  
    - Otherwise, if you do not have Intel Parallel Studio, open an instance of the **command-line interface** that comes with Microsoft Visual Studio from the list of programs in Windows start menu. This is simply a Windows command prompt that has all the necessary compiler variables and paths predefined in it.  

- **Build the ParaMonte example**:  
    Buid the example via Intel Parallel Studio command-line interface,  
    ```
    build.bat  
    ```
    The build script will automatically detect whether a parallel application has be built. By default, the name of the output executable is `runExample.exe`.  

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

### Linux  

**Environment**  


**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/pm) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
