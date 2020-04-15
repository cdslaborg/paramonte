%**********************************************************************************************************************************
%**********************************************************************************************************************************
%
%  ParaMonte: plain powerful parallel Monte Carlo library.
%
%  Copyright (C) 2012-present, The Computational Data Science Lab
%
%  This file is part of ParaMonte library. 
%
%  ParaMonte is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, version 3 of the License.
%
%  ParaMonte is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
%
%**********************************************************************************************************************************
%**********************************************************************************************************************************
classdef paramonte < handle
%   This is the Python interface to ParaMonte: Plain Powerful Parallel Monte Carlo library.
%   
%   What is ParaMonte?
%   ==================
%   
%   ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical 
%   objective functions of arbitrary-dimensions, in particular, the posterior distributions 
%   of Bayesian models in data science, Machine Learning, and scientific inference, with the 
%   design goal of unifying the 
%   
%       **automation** (of Monte Carlo simulations), 
%       **user-friendliness** (of the library), 
%       **accessibility** (from multiple programming environments), 
%       **high-performance** (at runtime), and 
%       **scalability** (across many parallel processors).  
%   
%   For more information on the installation, usage, and examples, visit: 
%   
%       https://www.cdslab.org/paramonte  
%   
%   The routines currently supported by the Python interface of ParaMonte include:
%   
%       ParaDRAM
%       ========
%           
%           Parallel Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo Sampler.
%   
%           EXAMPLE SERIAL USAGE
%           ====================
%   
%           Copy and paste the following code enclosed between the 
%           two comment lines in your python/ipython/jupyter session
%           (make sure the indentation of the pasted lines is correct):
%   
%   ##################################
%   import paramonte as pm
%   import numpy as np
%   def getLogFunc(Point):
%       # return the log of the standard multivariate 
%       # Normal density function with ndim dimensions
%       return -0.5 * np.sum( np.double( Point )**2 )
%   pmpd = pm.ParaDRAM()
%   pmpd.runSampler ( ndim = 3
%                   , getLogFunc = getLogFunc
%                   )
%   ##################################
%   
%           EXAMPLE PARALLEL USAGE
%           ======================
%   
%           Copy and paste the following code enclosed between the 
%           two comment lines in your python/ipython/jupyter session
%           (make sure the indentation of the pasted lines is correct):
%   
%   ##################################
%   with open("main.py", "w") as file:
%       file.write  ('''
%   import paramonte as pm
%   import numpy as np
%   def getLogFunc(Point):
%       # return the log of the standard multivariate 
%       # Normal density function with ndim dimensions
%       return -0.5 * np.sum( np.double( Point )**2 )
%   pmpd = pm.ParaDRAM()
%   pmpd.runSampler ( ndim = 3
%                   , getLogFunc = getLogFunc
%                   , mpiEnabled = True
%                   )
%   ''')
%   ##################################
%   
%           This will generate a main.py Python script file in the current
%           working directory of your Python session. Now, you can execute 
%           this Python script file (main.py) in parallel in two ways:
%           
%               1.  from inside ipython or jupyter: type the following,
%   
%                      !mpiexec -n 3 python main.py
%   
%               2.  outside of Python environment, 
%                   from within a Bash shell (on Linux or Mac) or,
%                   from within an Anaconda command prompt on Windows,
%                   type the following,
%   
%                      mpiexec -n 3 python main.py
%   
%               NOTE: On Windows platform, if you are using Intel MPI library,
%               NOTE: you may also specify the extra flag -localonly to run only
%               NOTE: on one node, but also to avoid the use of Hydra service and
%               NOTE: its registration. If you are not on a Windows cluster, (e.g., 
%               NOTE: you are using your personal device), then we recommend 
%               NOTE: specifying this flag.
%   
%               In both cases above, the script 'main.py' will run on 3 processors.
%               Feel free to change the number of processors to any number desired.
%               But do not request more than the number of physical cores on your system.
%
    properties (Access = public)
        authors             = []
        credits             = []
        version             = []
        rootPath            = []
        ParaDRAM            = []
    end
%
%    properties (Access = public, Hidden)
%        spec                = []
%        Stats               = ParaDRAM_Statistics_class()
%        Chain               = []
%        RefinedChain        = []
%    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = paramonte(varargin)

            % check interface type

            errorOccurred = false;
            matlabKernelEnabled = false;
            if nargin==1
                if isa(varargin{1},'char')
                    if strcmp(lower(varargin{1}),'matlab')
                        matlabKernelEnabled = true;
                    else
                        errorOccurred = true;
                    end
                elseif isa(varargin{1},'string')
                    if strcmp(lower(varargin{1}),"matlab")
                        matlabKernelEnabled = true;
                    else
                        errorOccurred = true;
                    end
                end
            elseif nargin~=0
                errorOccurred = true;
            end
            if errorOccurred
                disp( [ newline, 'The paramonte class constructor takes at most one argument of value ''matlab''. you have entered:', newline ] )
                disp(varargin)
                disp(   [ newline, 'Pass the input value ''matlab'' only if you know what it means. Otherwise, do not pass any input values. ' ...
                        , 'ParaMonte will properly set things up for you.', newline ...
                        ])
                disp('');
                error(['Gracefully exiting the ParaMonte constructor.']);
            end

            self.rootPath = mfilename('fullpath'); [self.rootPath,~,~] = fileparts(self.rootPath);
            addpath(genpath(self.rootPath));

            self.authors = "The Computational Data Science Lab @ The University of Texas";
            self.credits = "Peter O'Donnell Fellowship";

            versionFilePath = fullfile(self.rootPath,'.VERSION');
            versionFileID = fopen(versionFilePath);
            self.version = fgetl(versionFileID);

            if matlabKernelEnabled
                self.ParaDRAM = ParaDRAM_class;
            else
                self.ParaDRAM = ParaDRAM;
            end

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function result = getVersion(self)
            result = ["ParaMonte Python Interface Version " + self.version];
        end


    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end % methods (dynamic)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Static)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function result = helpme()
            result = help(paramonte);
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods (Static)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end