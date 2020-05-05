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

The instructions for building the MATLAB interface to the ParaMonte library are identical to the C interface, except that, 

    - the library type at the time of build must be set to `dynamic` (as opposed to `static`),
    - the exception handling is also enabled for MATLAB serial applications to avoid unexpected complete shut-downs of MATLAB.

## Usage  

**Environment**  

-   **Install MATLAB distribution >2018**: The MATLAB interface to the ParaMonte library minimally requires a recent MATLAB installation on your system. If you intend to build the MATLAB interface to ParaMonte library from scratch on your system, then you also need to install the MATLAB compilers component on your system as well.  

-   **Optionally install the MATLAB runtime libraries >2018**: If you intend to compile your MATLAB applications to run with ParaMonte, you will have to also download and install the MATLAB runtime libraries. Since the only way to run MATLAB in parallel is to compile your MATLAB files to generate executables, the MATLAB runtime libraries are required for parallel ParaMonte simulations.  

-   **Optionally install a compatible MPI library**: For parallel simulations (via MPI), you will also need an MPI library already installed on your system. We recommend Intel MPI library which is available for free from their website. On macOS, OpenMPI can be used as Intel MPI library currently (as of January 2020) does not support macOS.  

**Running the ParaMonte example**  

-   Open a MATLAB interface.  

-   For **serial** simulations, simply type the name of the MATLAB main script, which, by default, should be `main.m`.  
-   For **parallel** applications, first make sure that the input property `mpiEnabled` of your ParaMonte sampler is set to `true`, for example, `pmpd.mpiEnabled = true;` before calling the `runSampler()` method of your ParaMonte sampler object. Then call the main_mpi MATLAB file on the MATLAB command line,  

**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/pm) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
