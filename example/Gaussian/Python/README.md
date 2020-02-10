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

The instructions for building the Python interface of ParaMonte library are identical to the C interface, except that, the library type at the time of build must be set to `dynamic` (as opposed to `static`).

## Usage  

**Environment**  

- **Install Anaconda distribution of Python 3**: The Python interface of ParaMonte minimally requires the following Python libraries: `numpy`, `scipy`, `pandas`, `seaborn`. The Anaconda distribution of Python automatically comes with all of these packages when installed on your system.  

- **Optionally install a compatible MPI library**: For parallel simulations (via MPI), you will also need an MPI library already installed on your system. We recommend Intel MPI library which is avaialable for free from their website. On macOS, OpenMPI can be used as Intel MPI library currently (as of January 2020) does not support macOS.  

**Running the ParaMonte example**  

- Open an Anaconda command-line interface, or `jupyter` notebook.  

- For **serial** simulations, simply type the name of the Python main script, which, by default, should be `main.py`.  
- For **parallel** applications, first make sure that the input parameter `mpiEnabled = True` is passed to the `runSampler()` method of your ParaMonte sampler object. Then call the main Python file on the command line with the following syntax,  
    ```
    mpiexec -np NUM_PROCESSES runExample.exe
    ```
    where `NUM_PROCESSES` represents the number of processes on which the application will run.  

**For more information**, visit [cdslab.org/pm](https://www.cdslab.org/pm) or contact Amir Shahmoradi: [shahmoradi@utexas.edu](mailto:"shahmoradi@utexas.edu")  
