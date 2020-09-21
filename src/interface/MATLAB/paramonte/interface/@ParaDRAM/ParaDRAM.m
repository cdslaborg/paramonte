%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaDRAM - This is the ParaDRAM class for generating instances of serial and parallel
%   Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo
%   sampler of the ParaMonte library.
%   
%   Once you set the desired attributes to the desired values via the `spec`
%   component call the ParaDRAM sampler via the object's method runSampler().
%
%   Parameters
%   ----------
%
%       platform
%
%           A MATLAB struct containing the platform logical values (isMacOS, isWin32, isLinux).
%           Note that ll class attributes can be set after an instance
%           is returned by the constructor.
%
%   Attributes
%   ----------
%
%       See below for information on the attributes (properties).  
%
%   Methods
%   -------
%
%       See below for information on the methods.  
%   
%   Returns
%   -------
%
%       Object of class ParaDRAM
%
%   Minimal serial example
%   ----------------------
%
%       Here is a MATLAB script main.m for a serial ParaDRAM simulation.
%       Copy and paste this code into your MATLAB command line to run it:
%
%       Copy and paste the following code enclosed between the
%       two comment lines in your MATLAB session:
%
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           pmlibRootDir = './'; % if needed, change this path to the ParaMonte library root directory
%           addpath(genpath(pmlibRootDir));
%           pm = paramonte();
%           pmpd = pm.ParaDRAM();
%           pmpd.runSampler ( 4                 ... number of dimensions of the objective function
%                           , @(x) -sum(x.^2)   ... the natural log of the objective function
%                           );
%           pmpd.readChain();
%           pmpd.chainList{1}.plot.grid.plot();
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       The mathematical objective function in the above example is a
%       is a multivariate Normal distribution centered at the origin,
%       whose natural logarithm is returned by the lambda (Anonymous)
%       function defined as a function handle input to the ParaDRAM
%       sampler.
%
%       Running this code will generate a set of simulation output files (in the current
%       working directory of MATLAB). Among these, the file suffixed with "_report.txt" 
%       contains the full description of all input specifications of the ParaDRAM 
%       simulation as well as other information about the simulation results.
%
%   Parallel simulations
%   --------------------
%
%       0.  Ensure you need and will get a speedup by running the ParaDRAM sampler in parallel.
%           Typically, if a single evaluation of the objective function takes much longer than
%           a few milliseconds, your simulation may then benefit from the parallel simulation.
%
%       1.  First, ensure the required MPI libraries are installed on your System:
%           (You can skip this step if you know that you already have 
%           a compatible MPI library installed on your system). 
%           On the MATLAB command line type the following, 
%
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               pmlibRootDir = './'; % if needed, change this path to the ParaMonte library root directory
%               addpath(genpath(pmlibRootDir));
%               pm = paramonte();
%               pm.verify();
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           This will verify the existence of a valid MPI library on your system and,
%           if missing, will install the MPI library on your system (with your permission).
%
%       2.  Once the MPI installation is verified, 
%           copy and paste the following code enclosed 
%           between the two comment lines in your MATLAB session:
%
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               fid = fopen("main_mpi.m", "w");
%               sourceCode = ...
%               "pmlibRootDir = './'; % if needed, change this path to the ParaMonte library root directory" + newline + ...
%               "addpath(genpath(pmlibRootDir));" + newline + ...
%               "pm = paramonte();" + newline + ...
%               "pmpd = pm.ParaDRAM();" + newline + ...
%               "pmpd.mpiEnabled = true;" + newline + ...
%               "pmpd.runSampler ( 4                 ... number of dimensions of the objective function" + newline + ...
%               "                , @(x) -sum(x.^2)   ... the natural log of the objective function" + newline + ...
%               "                );";
%               fprintf( fid, "%s\n", sourceCode );
%               fclose(fid);
%               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       3.  This will generate a main_mpi.m MATLAB script file in the current
%           working directory of your MATLAB session. Now, you can execute
%           this MATLAB script file (main_mpi.m) in parallel. To do so, 
%           you need to call MATLAB on a command-line, out of MATLAB.
%
%               a.  On Windows: 
%
%                   From within command prompt that recognizes both MATLAB and mpiexec, 
%                   (ideally, the Intel Parallel Studio's command-prompt)
%                   type the following,
%
%                       mpiexec -localonly -n 3 matlab -batch "main_mpi"
%
%                   NOTE
%
%                       In the above MPI launcher commands for Windows OS, 
%                       we assumed that you would be using the Intel MPI library, hence, 
%                       the reason for the extra flag -localonly. This flag runs the parallel 
%                       code only on one node, but in doing so, it avoids the use of Hydra service 
%                       and its registration. If you are not on a Windows cluster, (e.g., you are 
%                       using your personal device), then we recommend specifying this flag.
%
%
%               b.  On macOS/Linux: 
%
%                   From within a Bash terminal that recognizes both MATLAB and mpiexec, 
%                   type the following,
%
%                       mpiexec -n 3 matlab -batch "main_mpi"
%
%           NOTE
%
%               In both cases in the above, the script 'main_mpi.m' will run on 3 processors.
%               Feel free to change the number of processors to any number desired. But do not 
%               request more than the available number of physical cores on your system.
%
%       WARNING
%
%           Do not add postprocessing codes (such as reading and plotting the output samples)
%           in your parallel code. There is no point in doing so, since MATLAB will run in `-batch`
%           mode for parallel simulations, disabling all plotting capabilities. Moreover, if you read
%           and postprocess the output files in parallel mode, the task will be done by all of the parallel
%           processes, potentially overwriting each others external activities.
%           Only perform the sampling (by calling the sampler routine) in parallel mode.
%
%   ParaDRAM Simulation Attributes
%   ------------------------------
%
%       The ParaDRAM simulation specifications have lengthy comprehensive descriptions
%       that appear in full in the output report files of every ParaDRAM simulation.
%
%       The best way to learn about individual ParaDRAM simulation attributes
%       is to a run a minimal serial simulation as given in the above.
%       You can also use the helpme() method of the spec component:
%
%           pm = paramonte();
%           pmpd = pm.ParaDRAM();
%           pmpd.spec.helpme();
%
%       See also: https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef ParaDRAM < ParaMonteSampler

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = ParaDRAM(platform,website)
            self = self@ParaMonteSampler(platform,website);
            self.spec = SpecDRAM(self.website.home.url);
            self.methodName = "ParaDRAM";
            self.method.isParaDRAM = true;
            self.Err.prefix = self.methodName;
            self.Err.resetEnabled = false;
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        runSampler(self,ndim,getLogFunc,varargin)
        [markovChainList] = readMarkovChain(self,varargin)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef ParaDRAM < ParaMonteSampler