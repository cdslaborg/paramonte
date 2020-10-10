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
%   paramonte - This is the MATLAB interface to ParaMonte: Plain Powerful Parallel Monte Carlo library.
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
%   To get quick help on the paramonte class, navigate to
%   the ParaMonte root directory and, type the following commands
%   enclosed between the two comment lines in your MATLAB session,
%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       pmlibRootDir = './'; % if needed, change this path to the ParaMonte library root directory
%       addpath(genpath(pmlibRootDir));
%       pm = paramonte();
%       doc pm
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Parameters
%   ----------
%
%       kernelType (optional)
%
%           An optional string with only one possible value "matdram".
%           If specified, the MATLAB implementation of the ParaMonte
%           kernel routines will be used. Currently, only the ParaDRAM
%           routine has an equivalent implementation in pure MATLAB.
%           Keep in mind that the kernel routines implemented in MATLAB
%           are typically 10 times slower than the equivalent routines
%           implemented in Fortran.
%           If not provided, the default Fortran kernel routines will be used.
%
%   Attributes
%   ----------
%
%       See below for information on the attributes (properties).
%
%   Methods
%   -------
%
%       See below for information on the methods and sampler constructors.
%
%   Naming conventions
%   ------------------
%
%       The camelCase naming style is used throughout the entire ParaMonte library, across
%       all programming languages: C/Fortran/Julia/MATLAB/Python
%
%       All simulation specifications start with a lowercase letter, including
%       scalar/vector/matrix int, float, string, or boolean variables.
%
%       The name of any variable that represents a vector of values is suffixed with "Vec",
%       for example: startPointVec, domainLowerLimitVec, ...
%
%       The name of any variable that represents a matrix of values is suffixed with "Mat",
%       for example: proposalStartCorMat, ...
%
%       The name of any variable that represents a list of varying-size values is suffixed
%       with "List", for example: variableNameList, ...
%
%       All functions or class methods begin with a lowercase verb.
%
%       Significant attempt has been made to end all boolean variables with a passive verb,
%       such that the full variable name virtually forms an English-language statement
%       that should be either True or False, set by the user.
%
%   Tips
%   ----
%
%       When running ParaMonte samplers, in particular on multiple cores in parallel,
%       it would be best to close any such aggressive software/applications as
%       Dropbox, ZoneAlarm, ... that can interfere with your ParaMonte
%       simulation output files, potentially causing the sampler to
%       crash before successful completion of the simulation.
%       These situations should however happen only scarcely.
%
%       On Windows systems, when restarting an old interrupted ParaDRAM simulation,
%       ensure your MATLAB session is also restarted before the simulation restart.
%       This may be needed as Windows sometimes locks access to some or all of the
%       simulation output files.
%
%       To unset an already-set input simulation specification, simply set the
%       simulation attribute to empty double `[]` or re-instantiate the object.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef paramonte %< dynamicprops

    properties (Access = public)
        %
        %       website
        %           A structure containing some web addresses relevant to the ParaMonte library
        website         = [];
        %
        %       authors
        %           A string containing the ParaMonte library authors
        authors         = [];
        %
        %       credits
        %           A string containing the the acknowledgment statement
        credits         = [];
        %
        %       version
        %           A structure containing the ParaMonte library version information
        %
        %           interface
        %               An object of class Version containing the
        %               ParaMonte library interface version get() and dump() methods.
        %               Example usage:
        %
        %                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   pm = paramonte();
        %                   pm.version.interface.get()  % prints the full interface version info
        %                   pm.version.interface.dump() % prints the only interface version
        %                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %           kernel
        %               An object of class Version containing the
        %               ParaMonte library kernel version get() and dump() methods.
        %               Example usage:
        %
        %                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   pm = paramonte();
        %                   pm.version.kernel.get()     % prints the full kernel version info
        %                   pm.version.kernel.dump()    % prints the only kernel version
        %                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        version         = [];
        %
        %           Parallel Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo Sampler.
        %
        %           To see the description and an example usage of the ParaDRAM routine,
        %           type the following commands enclosed between the two comment lines
        %           in your MATLAB session:
        %
        %               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %               pmlibRootDir = './'; % if needed, change this path to the ParaMonte library root directory
        %               addpath(genpath(pmlibRootDir));
        %               pm = paramonte();
        %               pmpd = pm.ParaDRAM();
        %               doc pmpd
        %               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ParaDRAM        = [];
    end

    properties(Hidden)
        %
        %           A MATLAB struct with three logical components:
        %
        %               isWin32 : true if the operating system is Windows
        %               isLinux : true if the operating system is Linux
        %               isMacOS : true if the operating system is macOS (Darwin)
        %
        platform = struct();
    end

    properties(Access = protected, Hidden)
        Err = Err_class();
        bashrcContentsBeforeParaMonteInstall;
        verificationStatusFilePath;
        buildInstructionNote;
        matdramKernelEnabled;
        pmInstallFailed;
        objectName;
        prereqs;
        isGUI;
        names;
        path;
    end

    properties (Access = protected)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = paramonte(varargin)

            self.path = struct();
            self.path.home = getHomePath();
            self.path.root = mfilename('fullpath');
            [self.path.root,~,~] = fileparts(self.path.root);
            self.path.root = string( getFullPath(fullfile(self.path.root,"..",".."),'lean') );
            addpath(genpath(self.path.root),'-begin');
            self.path.lib = string(fullfile(self.path.root, "lib"));
            self.path.auxil = string(fullfile(self.path.root, "auxil"));

            for versionType = ["interface","kernel"]
                self.version.(versionType) = Version_class(self.path.auxil,versionType);
            end

            %self.prereqs = struct();

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% set up website
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.website = struct();
            self.website.home = struct();
            self.website.home.url = "https://www.cdslab.org/paramonte";
            self.website.home.install = struct();
            self_website_home_install_url = self.website.home.url + "/notes/installation";

            %%%% installation Linux

            self.website.home.overview = struct();
            self_website_home_overview_url = self.website.home.url + "/notes/overview";
            self.website.home.overview.preface = struct();
            self.website.home.overview.changes = struct();
            self.website.home.overview.preface.url = self_website_home_overview_url + "/preface";
            self.website.home.overview.changes.kernel = struct();
            self.website.home.overview.changes.python = struct();
            self.website.home.overview.changes.matlab = struct();
            self.website.home.overview.changes.kernel.url = self_website_home_overview_url + "/paramonte-kernel-release-notes";
            self.website.home.overview.changes.python.url = self_website_home_overview_url + "/paramonte-python-release-notes";
            self.website.home.overview.changes.matlab.url = self_website_home_overview_url + "/paramonte-matlab-release-notes";

            %%%% installation Linux

            self.website.home.install.linux = struct();
            self.website.home.install.linux.url = self_website_home_install_url + "/linux";

            %%%% installation Windows

            self.website.home.install.windows = struct();
            self.website.home.install.windows.url = self_website_home_install_url + "/windows";

            %%%% installation MATLAB

            self.website.home.install.matlab = struct();
            self.website.home.install.matlab.url = self_website_home_install_url + "/matlab";

            %%%% installation Python

            self.website.home.install.python = struct();
            self.website.home.install.python.url = self_website_home_install_url + "/python";

            %%%% installation macOS

            self.website.home.install.macos = struct();
            self.website.home.install.macos.url = self_website_home_install_url + "/macos";
            self.website.home.install.macos.prereqs = struct();
            self.website.home.install.macos.prereqs.url = self.website.home.install.macos.url + "/#the-compile-time-and-runtime-prerequisites";
            self.website.home.install.macos.prereqs.cmd = struct();
            self.website.home.install.macos.prereqs.cmd.url = self.website.home.install.macos.url + "/#prereqs-install";

            %%%% MATLAB examples

            self.website.home.examples = struct();
            self_website_home_examples_url = self.website.home.url + "/notes/examples";
            self.website.home.examples.matlab = struct();
            self.website.home.examples.matlab.jupyter = struct();
            self.website.home.examples.matlab.postprocess = struct();
            self.website.home.examples.matlab.jupyter.url = self_website_home_examples_url + "/matlab/jupyter";
            self.website.home.examples.matlab.postprocess.url = self_website_home_examples_url + "/matlab/postprocess";

            %%%% Python examples

            self.website.home.examples = struct();
            self_website_home_examples_url = self.website.home.url + "/notes/examples";
            self.website.home.examples.python = struct();
            self.website.home.examples.python.jupyter = struct();
            self.website.home.examples.python.postprocess = struct();
            self.website.home.examples.python.jupyter.url = self_website_home_examples_url + "/python/jupyter";
            self.website.home.examples.python.postprocess.url = self_website_home_examples_url + "/python/postprocess";

            %%%% Python API

            self.website.home.api = struct();
            self_website_home_api_url = self.website.home.url + "/notes/api";
            self.website.home.api.python = struct();
            self.website.home.api.python.url = self_website_home_api_url + "/python/autoapi/paramonte";

            %%%% ParaDRAM

            self.website.home.usage = struct();
            self_website_home_usage_url = self.website.home.url + "/notes/usage";
            self.website.home.usage.paradram = struct();
            self_website_home_usage_paradram_url = self_website_home_usage_url + "/paradram";
            self.website.home.usage.paradram.quickstart = struct();
            self.website.home.usage.paradram.quickstart.url = self_website_home_usage_paradram_url + "/interface";
            self.website.home.usage.paradram.input = struct();
            self.website.home.usage.paradram.input.url = self_website_home_usage_paradram_url + "/input";
            self.website.home.usage.paradram.specifications = struct();
            self.website.home.usage.paradram.specifications.url = self_website_home_usage_paradram_url + "/specifications";
            self.website.home.usage.paradram.restart = struct();
            self.website.home.usage.paradram.restart.url = self_website_home_usage_paradram_url + "/restart";
            self.website.home.usage.paradram.output = struct();
            self.website.home.usage.paradram.output.url = self_website_home_usage_paradram_url + "/output";

            %%%% GitHub issues

            self.website.github = struct();
            self.website.github.url = "https://github.com/cdslaborg/paramonte";
            self.website.github.issues = struct();
            self.website.github.issues.url = "https://github.com/cdslaborg/paramonte/issues";
            self.website.github.release = struct();
            self.website.github.release.url = self.website.github.url + "/releases";
            self.website.github.release.latest = struct();
            self.website.github.release.latest.url = self.website.github.release.url + "/latest";
            self.website.github.archive = struct();
            self_website_github_archive_url = self.website.github.url + "/archive";
            self.website.github.archive.master = struct();
            self.website.github.archive.master.zip = struct();
            self.website.github.archive.master.tar = struct();
            self.website.github.archive.master.zip.url = self_website_github_archive_url + "/master.zip";
            self.website.github.archive.master.tar.url = self_website_github_archive_url + "/master.tar.gz";

            %%%% GitHub examples

            self.website.github.examples = struct();
            self.website.github.examples.url = "https://github.com/cdslaborg/paramontex";

            %%%% Intel MPI

            self.website.intel = struct();
            self.website.intel.mpi = struct();
            self.website.intel.mpi.home = struct();
            self.website.intel.mpi.home.url = "https://software.intel.com/en-us/mpi-library";

            %%%% Intel MPI Windows

            self.website.intel.mpi.windows = struct();
            self.website.intel.mpi.windows.url = "https://software.intel.com/en-us/get-started-with-mpi-for-windows";

            %%%% OpenMPI

            self.website.openmpi = struct();
            self.website.openmpi.home = struct();
            self.website.openmpi.home.url = "https://www.open-mpi.org/";

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            self.isGUI = isGUI(); 

            self.platform.isWin32 = ispc;
            self.platform.isMacOS = ismac;
            self.platform.isLinux = isunix && ~ismac;
            if self.platform.isWin32; self.platform.osname = "Windows"; end
            if self.platform.isMacOS; self.platform.osname = "Darwin"; end % do NOT set this to macOS !
            if self.platform.isLinux; self.platform.osname = "Linux"; end

            self.platform.systemInfoFilePrefix = fullfile(self.path.auxil, ".systemInfo_");
            self.platform.systemInfoFilePath = self.platform.systemInfoFilePrefix + string(strrep(date,"-","_"));
            if ~isfile(self.platform.systemInfoFilePath)
                delete(self.platform.systemInfoFilePrefix + "*");
                if self.platform.isWin32
                    cmd = "systeminfo";
                end
                if self.platform.isMacOS
                    cmd = "uname -a; sysctl -a | grep machdep.cpu";
                end
                if self.platform.isLinux
                    cmd = "uname -a; lscpu";
                end
                [errorOccurred, self.platform.sysinfo] = system(cmd);
                if errorOccurred
                    warning ( newline ...
                            + "The ParaMonte library failed to fetch the system information on your system. skipping..." ...
                            + newline ...
                            );
                end
                fid = fopen(self.platform.systemInfoFilePath,'wt');
                fprintf(fid, "%s", self.platform.sysinfo);
                fclose(fid);
            end

            self.authors = "The Computational Data Science Lab @ The University of Texas";
            self.credits = "Peter O'Donnell Fellowship / Texas Advanced Computing Center";

            self.verificationStatusFilePath = fullfile( self.path.auxil, ".verificationEnabled" );
            self.buildInstructionNote   = "If your platform is non-Windows and is compatible with GNU Compiler Collection (GCC)," + newline ...
                                        + "you can also build the required ParaMonte kernel's shared object files on your system" + newline ...
                                        + "by calling ParaMonte module's build() function from within your MATLAB environment, like:" + newline + newline ...
                                        + "    pm = paramonte();" + newline ...
                                        + "    pm.build()";

            self.names.paramonte = "ParaMonte";
            self.names.paradram = "ParaDRAM";
            self.names.paranest = "ParaNest";
            self.names.paratemp = "ParaTemp";

            self.Err.marginTop = 1;
            self.Err.marginBot = 1;
            self.Err.prefix = self.names.paramonte;
            self.Err.resetEnabled = false;

            % check interface type

            errorOccurred = false;
            matlabKernelName = "matdram";
            matlabInterfName = "interface";
            if nargin==0
                self.matdramKernelEnabled = true; if isfolder(self.path.lib); self.matdramKernelEnabled = false; end
            elseif nargin==1
                if isa(varargin{1},"char") || isa(varargin{1},"string")
                    if strcmpi(string(varargin{1}),matlabKernelName)
                        self.matdramKernelEnabled = true;
                    elseif strcmpi(string(varargin{1}),matlabInterfName)
                        self.matdramKernelEnabled = false;
                    else
                        errorOccurred = true;
                    end
                end
            else
                errorOccurred = true;
            end
            if errorOccurred
                self.Err.msg    = "The paramonte class constructor takes at most one argument of value """ ...
                                + matlabKernelName + """ or """ + matlabInterfName + """. You have entered:" + newline + newline ...
                                ... + "    " + string(strrep(join(string(varargin)," "),'\','\\')) ...
                                + "    " + join(string(varargin)," ") ...
                                + "Pass the input value """ + matlabKernelName + """ only if you know what it means. " ...
                                + "Otherwise, do not pass any input values. ParaMonte will properly set things up for you.";
                self.Err.abort()
            end

            % set ParaDRAM sampler

            if self.matdramKernelEnabled
                self.ParaDRAM = ParaDRAM_class(self.platform,self.website,self.version.interface.dump);
            else
                self.ParaDRAM = ParaDRAM(self.platform,self.website);
            end

            % verify prereqs

            self.verify("reset",false);

            if self.matdramKernelEnabled; return; end

            if ~self.platform.isWin32
                LD_LIBRARY_PATH = getenv("LD_LIBRARY_PATH");
                if isempty(LD_LIBRARY_PATH); LD_LIBRARY_PATH = ""; end
                setenv("LD_LIBRARY_PATH",self.path.lib+":"+LD_LIBRARY_PATH);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function verify(self,varargin)
            %
            %   checks (or rechecks) the requirements of the installed ParaMonte library
            %
            %   Parameters
            %   ----------
            %
            %       reset
            %           boolean whose default value is true. If true,
            %           a thorough verification of the existence of the required
            %           libraries will performed, as if it is the first ParaMonte
            %           module import.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Example usage
            %   -------------
            %
            %       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       pm = paramonte()
            %       pm.verify()                 % default
            %       pm.verify(false)            % default, same as above
            %       pm.verify("reset",false)    % default, same as above
            %       pm.verify("reset",true)     % resets and performs all checks
            %       pm.verify(true)             % resets and performs all checks
            %       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            self.objectName = inputname(1);

            setenv('PATH', ['/usr/local/bin:',getenv('PATH')]);

            invalidArgumentProvided = false;
            if nargin==1
                reset = true;
            elseif nargin==2
                if isa(varargin{1},"logical")
                    reset = varargin{1};
                else
                    invalidArgumentProvided = true;
                end
            elseif nargin==3
                if ( isa(varargin{1},"string") || isa(varargin{1},"char") ) && ( strcmpi(varargin{1},"reset") && isa(varargin{2},"logical") )
                    reset = varargin{2};
                else
                    invalidArgumentProvided = true;
                end
            end
            if invalidArgumentProvided
                fullMethodName = self.objectName + ".verify";
                error   ( newline ...
                        + "Invalid usage of " + fullMethodName + "() detected. Possible valid usages are:" + newline + newline ...
                        + "    " + fullMethodName + "()" + newline ...
                        + "    " + fullMethodName + "(true) % : the verification status in the ParaMonte cache will be reset to the default initial value." + newline ...
                        + "    " + fullMethodName + "('reset',true) % same as above." ...
                        );
            end

            % require MATLAB >2017a

            matlabVersion = version("-release");
            matlabVersion = str2double(matlabVersion(1:4));
            if matlabVersion < 2017
                error   ( newline ...
                        + "MATLAB R" + string(matlabVersion) + "a or newer is required for proper ParaMonte functionality. " ...
                        + "Please install a compatible version of MATLAB." ...
                        + newline ...
                        );
            end

            % ensure messages are printed only for the first-time import

            if reset; self.writeVerificationStatusFile("True"); end

            try
                fid = fopen(self.verificationStatusFilePath);
                verificationEnabledString = fgetl(fid);
                fclose(fid);
            catch
                self.Err.msg    = "The Integrity of the ParaMonte library appears to have been compromised. " + newline ...
                                + "The following file is missing or not readable in the ParaMonte library: " + newline + newline ...
                                + "    " + self.verificationStatusFilePath + newline + newline ...
                                + "To fix this issue, download a fresh clone of the ParaMonte MATLAB library from, " + newline + newline ...
                                + "    " + href(self.website.github.release.url) ...
                                ;
                self.Err.abort();
            end
            if verificationEnabledString=="False"
                verificationEnabled = false;
            elseif verificationEnabledString=="True" || verificationEnabledString=="Testing"
                verificationEnabled = true;
            else
                error   ( newline ...
                        + "The internal settings of the ParaMonte library appears to have been messed up " ...
                        + "potentially by the operating system, MATLAB, or some third party applications. " ...
                        + "Please reinstall ParaMonte by typing the following commands " ...
                        + "on a MATLAB-aware command-line interface:" + newline + newline ...
                        + "    pip uninstall paramonte" + newline ...
                        + "    pip install paramonte" ...
                        + newline ...
                        );
            end

            if verificationEnabled

                self.displayParaMonteBanner();

                % if MatDRAM, return

                if self.matdramKernelEnabled
                    self.Err.msg = "The ParaMonte MatDRAM library is ready to use!";
                    self.Err.note();
                    self.writeVerificationStatusFile("False");
                    return;
                end

                % ensure 64-bit architecture

                if strcmpi(getArch(),"x64") && (self.platform.isWin32 || self.platform.isLinux || self.platform.isMacOS)

                    % set up library path

                    if ~self.platform.isWin32
                        % save a copy of the original bashrc contents. will be used in installParaMonte().
                        self.bashrcContentsBeforeParaMonteInstall = getBashrcContents();
                        self.setupUnixPath();
                    end

                    % install library

                    self.installParaMonte();

                    % get the dependency list

                    self.prereqs = self.getPrereqs(false); % sets the self.prereqs component.

                    % search for the MPI library

                    mpi = self.findMPI();

                    if mpi.install.found && ~mpi.path.broken

                        self.writeVerificationStatusFile("False");

                    else

                        if mpi.install.found && mpi.path.broken
                            msg = "An MPI library installation appears to exist on your system, however, " ...
                                + "some components of the library appear to be missing, or the environmental path " ...
                                + "to the MPI library installation is corrupted. You can inspect the contents of the " ...
                                + "environmental path variable for potential path corruptions by typing," ...
                                + newline + newline ...
                                + "    " + "getenv(""PATH"")" ...
                                + newline + newline ...
                                + "on your MATLAB command line. If you or the ParaMonte library (on your behalf) " ...
                                + "have already successfully installed an MPI library on your system, " ...
                                + "you can safely ignore this warning and avoid further reinstallation of the MPI library. " ...
                                + "Otherwise, follow the instructions below to reinstall the MPI library. " ...
                                ;
                        else
                            msg = "The MPI runtime libraries for 64-bit architecture could not be detected on your system. " ...
                                + "The MPI runtime libraries are required for the parallel ParaMonte simulations. " ...
                                ;
                        end

                        self.Err.msg    = msg ...
                                        + "For Windows and Linux operating systems, you can download and install the Intel MPI runtime " ...
                                        + "libraries, free of charge, from Intel website, " ...
                                        + newline + newline ...
                                        + "    " + href(self.website.intel.mpi.home.url) ...
                                        + newline + newline ...
                                        + "For macOS (Darwin operating system), you can download and install the Open-MPI library. " ...
                                        + newline + newline ...
                                        + "    " + href(self.website.openmpi.home.url) ...
                                        + newline + newline ...
                                        + "Alternatively, the ParaMonte library can automatically install these library for you now. " ...
                                        + "If you don't know how to download and install the correct MPI runtime library version, " ...
                                        + "we strongly recommend that you let the ParaMonte library to install this library for you. " ...
                                        + "If so, ParaMonte will need access to the World-Wide-Web to download the library " ...
                                        + "and will need your administrative permission to install it on your system. " ...
                                        + "Therefore, if you have any active firewall on your system such as ZoneAlarm, " ...
                                        + "please make sure your firewall allows ParaMonte to access the Internet.";
                        self.Err.warn();

                        if verificationEnabledString=="True"

                            isYes = getUserResponse ( newline ...
                                                    + "    Do you wish to download and install the MPI runtime library" + newline ...
                                                    + "    (only needed for parallel simulations) on your system now (y/n)? " ...
                                                    );
                            if isYes
                                succeeded = self.installMPI();
                                if succeeded
                                    self.writeVerificationStatusFile("Testing");
                                else
                                    self.writeVerificationStatusFile("False");
                                end
                            else
                                self.Err.msg    = "Skipping the MPI library installation... " + newline ...
                                                + "It is now the user's responsibility to ensure an MPI runtime library" ...
                                                + "exists on the system for parallel simulations. " ...
                                             ...+ "If you ever wish to install MPI libraries via ParaMonte again, " ...
                                             ...+ "please try:" + newline + newline ...
                                             ...+ "    pm = paramonte();" + newline ...
                                             ...+ "    pm.verify();" + newline + newline ...
                                                + "For more information visit:" + newline + newline ...
                                                + "    " + href(self.website.home.url);
                                self.Err.note();
                                self.writeVerificationStatusFile("False");
                            end

                        end

                    end

                else

                    self.warnForUnsupportedPlatform();
                    succeeded = self.build();
                    if succeeded
                        self.writeVerificationStatusFile("Testing");
                    else
                        self.writeVerificationStatusFile("False");
                    end

                end

                self.Err.prefix = self.names.paramonte;
                self.Err.marginTop = 1;
                self.Err.marginBot = 1;
                self.Err.msg    = "To check for the ParaMonte or the MPI library installations status or, " + newline ...
                                + "to display the above messages in the future, use the following " + newline ...
                                + "commands on the MATLAB command-line: " + newline ...
                                + newline ...
                                + "    pm = paramonte();" + newline ...
                                + "    pm.verify();" ...
                                ;
                self.Err.note();

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function cite(self)
            disp(href(self.website.home.overview.preface.url + "/#how-to-acknowledge-the-use-of-the-paramonte-library-in-your-work"));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function checkForUpdate(self)

            currentVersionString = self.version.interface.dump();
            versionFileLineList = strsplit(webread("https://raw.githubusercontent.com/cdslaborg/paramonte/master/src/interface/MATLAB/.VERSION"), newline);
            latestVersionString = string(versionFileLineList{1});

            self.Err.prefix = self.names.paramonte;
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;

            while getVecLen(latestVersionString)
                try
                    latestVersionURL = "https://github.com/cdslaborg/paramonte/releases/tag/" + latestVersionString;
                    webread(latestVersionURL);
                    break;
                catch
                    latestVersionString = getPreviousVersion(latestVersionString);
                end
            end

            if getVecLen(latestVersionString)
                currentVersionTriplet = getVersionTriplet(currentVersionString);
                latestVersionTriplet = getVersionTriplet(latestVersionString);
                updateAvailable =   (latestVersionTriplet(1)  > currentVersionTriplet(1)) ...
                                ||  (latestVersionTriplet(1) == currentVersionTriplet(1) && latestVersionTriplet(2)  > currentVersionTriplet(2)) ...
                                ||  (latestVersionTriplet(1) == currentVersionTriplet(1) && latestVersionTriplet(2) == currentVersionTriplet(2) && latestVersionTriplet(3) > currentVersionTriplet(3));
            else
                updateAvailable = false;
            end

            if updateAvailable
                latestVersionFile = latestVersionURL + "/libparamonte_MATLAB";
                if self.platform.isWin32
                    latestVersionFile = latestVersionFile + ".zip";
                else
                    latestVersionFile = latestVersionFile + ".tar.gz";
                end
                self.Err.msg    = "A newer version (" + latestVersionString + ") of the ParaMonte library appears " + newline ...
                                + "to be available on the GitHub repository. The currently-installed version is: " + currentVersionString + newline ...
                                + "You can download the latest version of the ParaMonte MATLAB library from " + newline ...
                                + newline ...
                                + "    " + href(latestVersionFile) + newline ...
                                + newline ...
                                ;
            else
                self.Err.msg    = "It appears that you have the latest version of the ParaMonte library. " + newline;
            end

            self.Err.msg    = self.Err.msg ...
                            + "To see the latest changes to the ParaMonte MATLAB library, visit, " + newline ...
                            + newline ...
                            + "    " + href(self.website.home.overview.changes.matlab.url) ...
                            ;
            self.Err.note();

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (dynamic)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Access=public, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function installParaMonte(self)

            if isunix

                %pmLocalFileList = getFileNameList(self.path.lib);
                if self.platform.isLinux; filePath = "*so"; end
                if self.platform.isMacOS; filePath = "*dylib"; end
                fileList = dir(fullfile(self.path.lib,filePath));
                fileListLen = length(fileList);
                if fileListLen==0
                    self.Err.msg    = "Failed to locally detect the ParaMonte library files on your system. " ...
                                    + "The ParaMonte library folder appears to be empty. Please build a " ...
                                    + "fresh copy of the library or download it from, "  + newline + newline ...
                                    + "    " + href(self.website.github.release.url) ...
                                    ;
                    self.Err.abort();
                end

                self.pmInstallFailed = true;
                installRootDirList = ["/usr/local/lib64","/usr/lib64","/usr/local/lib","/usr/lib"];

                for installRootDir = installRootDirList
                    if isfolder(installRootDir)
                        errorOccurred = false;
                        %pmInstallDir = fullfile(installRootDir,"paramonte");
                        %if ~isfolder(pmInstallDir)
                        %    [status, errMsg, msgID] = mkdir(pmInstallDir);
                        %    if status~=1; errorOccurred = true; end
                        %end
                        %if ~errorOccurred
                            %filePath = "*";
                            dummy1 = self.Err.marginTop;
                            dummy2 = self.Err.marginBot;
                            self.Err.marginTop = 0;
                            self.Err.marginBot = 0;
                            for i = 1:fileListLen
                                fullFilePath = fullfile(fileList(i).folder,fileList(i).name);
                                %[status, errMsg] = system("sudo cp -rf " + fullFilePath + installRootDir + " &")
                                [status, errMsg, msgID] = copyfile(fullFilePath, installRootDir, "f");
                                if status==0 % && self.pmInstallFailed
                                    errorOccurred = true;
                                else
                                    self.pmInstallFailed = false;
                                    self.Err.msg = "The local installation of the ParaMonte library succeeded.";
                                    self.Err.note();
                                end
                            end
                            self.Err.marginTop = dummy1;
                            self.Err.marginBot = dummy2;
                            if errorOccurred && contains(installRootDir,"local")
                                self.Err.msg    = "An attempt to locally install the ParaMonte library on your system failed with the following message: " + newline  + newline ...
                                                + string(errMsg) + " Error flag: " + string(msgID) + newline  + newline ...
                                                + "Continuing at the risk of not being able to use the ParaMonte kernel samplers...";
                                self.Err.warn();
                            %else
                            %    self.pmInstallFailed = false;
                            %    break;
                            end
                            %[status, errMsg, msgID] = copyfile(fullfile(self.path.lib,filePath), installRootDir, "f");
                            %if status~=0; errorOccurred = true; end
                            %pmInstallFileList = getFileNameList(pmInstallDir);
                            %for localFile = pmLocalFileList
                            %    if ~any(strcmp(pmInstallFileList, localFile))
                            %        [status, errMsg, msgID] = copyfile(fullfile(self.path.lib,localFile), pmInstallDir, "f");
                            %        if status~=0; errorOccurred = true; break; end
                            %    end
                            %end
                        %end
                        % if errorOccurred && contains(installRootDir,"local")
                        %     self.Err.msg    = "An attempt to locally install the ParaMonte library on your system failed with the following message: " + newline  + newline ...
                        %                     + string(errMsg) + " Error flag: " + string(msgID) + newline  + newline ...
                        %                     + "Continuing at the risk of not being able to use the ParaMonte kernel samplers.";
                        %     self.Err.warn();
                        % %else
                        %     self.pmInstallFailed = false;
                        % %    break;
                        % end
                    end
                    if ~self.pmInstallFailed; break; end
                end

                self.Err.prefix = self.names.paramonte;
                self.Err.marginTop = 1;
                self.Err.marginBot = 1;
                if self.pmInstallFailed
                    self.Err.msg    = "Failed to locally install the ParaMonte library on this system. This is highly unusual. " ...
                                    + "If you have already successfully installed the ParaMonte library files locally on your system " ...
                                    + "through the following procedure, you can safely ignore this warning message and the subsequent " ...
                                    + "request for reinstallation given below. Otherwise, if you have administrator previlages on this system " ...
                                    + "(e.g., you are not using a supercomputer), please close and reopen MATLAB with ""sudo matlab"" " ...
                                    + "command from the Bash terminal. Then, navigate to the root directory of the ParaMonte library " ...
                                    + "to reinstall it via: " + newline ...
                                    + newline ...
                                    + "     pm = paramonte();" + newline ...
                                    + "     pm.verify();" + newline ...
                                    + newline ...
                                    + "This is essential to ensure the ParaMonte library files can be " ...
                                    + "written to the local installation disrectory. " ...
                                    + "If the local installation fails even after running MATLAB with administrator previlages, " ...
                                    + "please report this issue at: " ...
                                    + newline + newline ...
                                    + "    " + href(self.website.github.issues.url) ...
                                    ;
                    warningEnabled = true;
                    if self.platform.isLinux
                        if contains(self.bashrcContentsBeforeParaMonteInstall,self.path.lib)
                            warningEnabled = false;
                        else
                            self.Err.msg    = self.Err.msg ...
                                            + newline + newline ...
                                            ...+ "If you have already installed an older version of the ParaMonte library on this system in the past, " ...
                                            ...+ "and you used ""sudo matlab"" to install the library locally on your system, it is imperative that, " ...
                                            ...+ "You follow the same approach again to install the new version of the ParaMonte library. " ...
                                            + "Otherwise, as an alternative solution, you can follow the instructions given below," ...
                                            + newline + newline ...
                                            + self.getPathSetupMsg ...
                                            ;
                            question    = newline + "    Shall we quit this MATLAB session now so that you can " ...
                                        + newline + "    perform one of the suggested solutions in the above (y/n)? ";
                        end
                    else
                        question    = newline + "    Shall we quit this MATLAB session now so that you can " ...
                                    + newline + "    perform the suggested solution in the above (y/n)? "; % should be used only on macOS
                    end
                    if warningEnabled
                        self.Err.warn();
                        matlabExitRequested = getUserResponse( question );
                        if matlabExitRequested
                            exit;
                        else
                            self.Err.msg    = "Continuing at the risk of not being able to call the ParaMonte library kernel samplers...";
                            self.Err.warn();
                        end
                    end
                end

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function succeeded = build(self,varargin)

            succeeded = false;

            if nargin==1
                flags = "";
            elseif nargin==1
                flags = varargin{1};
            else
                error("invalid number of input arguments.");
            end

            self.Err.prefix = self.names.paramonte;
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;

            if self.platform.isWin32

                self.Err.msg    = "The ParaMonte library build on Windows Operating Systems (OS) requires " + newline ...
                                + "the installation of the following software on your system:" + newline + newline ...
                                + "    -- Microsoft Visual Studio (MSVS) (Community Edition >2017)" + newline ...
                                + "    -- Intel Parallel Studio >2018, which is built on top of MSVS" + newline + newline ...
                                + "If you don't have these software already installed on your system, " + newline ...
                                + "please visit the following page for the installation instructions:" + newline + newline ...
                                + "    " + href(self.website.home.url) + newline + newline ...
                                + "Follow the website's instructions to build the ParaMonte library on your system.";
                self.Err.note();
                return

            else

                self.Err.msg    = "You are requesting to build the ParaMonte kernel libraries on your system. " + newline ...
                                + "The kernel library build requires ParaMonte-compatible versions of the following " + newline ...
                                + "compilers and parallelism libraries to be installed on your system: " + newline + newline ...
                                + "    GNU compiler collection (GCC >8.3)" + newline ...
                                + "    MPI library (MPICH >3.2) on Linux OS or Open-MPI on Darwin OS" + newline ...
                                + "    OpenCoarrays >2.8" + newline + newline ...
                                + "The full installation of these software will require 4 to 5 Gb of free space " + newline ...
                                + "on your system (where the ParaMonte-MATLAB interface is already installed)." + newline ...
                                + "Note that the installation script is in Bash and therefore requires Bash shell." + newline ...
                                + "An existing recent installation of the GNU Compiler Collection (GCC) on your" + newline ...
                                + "system would be also highly desirable and will significantly cut the build time." + newline ...
                                + "Also, downloading the prerequisites requires access to the Internet." + newline ...
                                + "If you have an Internet firewall active on your system, please make sure to" + newline ...
                                + "turn it off before proceeding with the local installation of ParaMonte.";
                self.Err.note();

                buildEnabled = getUserResponse  ( newline + "    Do you wish to download and install the ParaMonte library" ...
                                                + newline + "    and its prerequisites on your system now (y/n)? " ...
                                                );

                if buildEnabled

                    if self.platform.isMacOS
                        self.buildParaMontePrereqsForMac();
                    end

                    pmGitTarPath = fullfile( self.path.lib, "master.tar.gz" );
                    pmGitTarPath = websave(pmGitTarPath,"https://github.com/cdslaborg/paramonte/archive/master.tar.gz");
                    pmGitRootDir = fullfile(self.path.lib, "paramonte-master");

                    try
                        untar(pmGitTarPath,self.path.lib);
                        pmGitInstallScriptPath = fullfile(pmGitRootDir, "install.sh");
                        if ~isfile(pmGitInstallScriptPath)
                            self.Err.msg    = "Internal error occurred." + newline ...
                                            + "Failed to detect the ParaMonte installation Bash script." + newline ...
                                            + "Please report this issue at " + newline + newline ...
                                            + "    " + href(self.website.github.issues.url) + newline + newline ...
                                            + "Visit, " ...
                                            + "    " + href(self.website.home.url) + newline + newline ...
                                            + "for instructions to build ParaMonte library on your system.";
                            self.Err.abort();
                        end
                    catch
                        self.Err.msg    = "Unzipping of the ParaMonte tarball failed." + newline ...
                                        + "Make sure you have tar software installed on your system and try again.";
                        self.Err.abort();
                    end

                    [errorOccurred, ~] = system("chmod +x " + string(strrep(pmGitInstallScriptPath,'\','\\')), "-echo");
                    if errorOccurred
                        self.Err.msg    = "The following action failed:" + newline + newline ...
                                        ... + "chmod +x " + string(strrep(pmGitInstallScriptPath,'\','\\')) + newline + newline ...
                                        + "chmod +x " + string(pmGitInstallScriptPath) + newline + newline ...
                                        + "skipping...";
                        self.Err.warn();
                    end

                    originalDir = cd(pmGitRootDir);

                    [errorOccurred1, ~] = system( ['find ', pmGitRootDir, ' -type f -iname \"*.sh\" -exec chmod +x {} \;'], "-echo" );
                    [errorOccurred2, ~] = system( pmGitInstallScriptPath + " --lang matlab --test_enabled true --exam_enabled false --yes-to-all " + flags, "-echo" );
                    if errorOccurred1 || errorOccurred2
                        self.Err.msg    = "The local installation of ParaMonte failed." + newline ...
                                        + "Please report this issue at " + newline + newline ...
                                        + "    " + href(self.website.github.issues.url);
                        self.Err.abort();
                    end

                    cd(originalDir);

                    % copy files to the library folder

                    matlabBinDir = string(fullfile( pmGitRootDir , "bin" , "MATLAB" , "paramonte" ));
                    fileList = dir( fullfile( matlabBinDir , "libparamonte_*" ) );

                    if isempty(fileList)

                        self.Err.msg    = "ParaMonte kernel libraries build and installation appears to have failed. " + newline ...
                                        + "You can check this path:" + newline + newline ...
                                        ... + string(strrep(matlabBinDir,'\','\\')) + newline + newline ...
                                        + matlabBinDir + newline + newline ...
                                        + "to find out if any shared objects with the prefix 'libparamonte_' have been generated or not." + newline ...
                                        + "Please report this issue at " + newline + newline ...
                                        + "    " + href(self.website.github.issues.url);
                        self.Err.abort();

                    else

                        self.Err.msg    = "ParaMonte kernel libraries build appears to have succeeded. " + newline ...
                                        + "copying the kernel files to the paramonte MATLAB module directory...";
                        self.Err.note();

                        for file = string({fileList(:).name})
                            %self.Err.msg = "file: " + string(strrep(file,'\','\\'));
                            self.Err.msg = "file: " + file;
                            self.Err.marginTop = 0;
                            self.Err.marginBot = 0;
                            self.Err.note();
                            shutil.copy(file, self.path.lib);
                        end

                        self.Err.msg = "ParaMonte kernel libraries should be now usable on your system.";
                        self.Err.marginTop = 1;
                        self.Err.marginBot = 1;
                        self.Err.note();

                        setupFilePath = fullfile( pmGitRootDir , "build", "prerequisites", "prerequisites", "installations", "opencoarrays", "2.8.0", "setup.sh" );

                        if isfile(setupFilePath)
                            bashrcContents = getBashrcContents();
                            setupFilePathCmd = "source " + setupFilePath;
                            if ~contains(bashrcContents,setupFilePathCmd)
                                [~, ~] = system( "chmod 777 ~/.bashrc", "-echo");
                                [~, ~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                                [~, ~] = system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte library local installation setup >>>' >> ~/.bashrc", "-echo" );
                                [~, ~] = system( "chmod 777 ~/.bashrc && echo '" + setupFilePathCmd + "' >>  ~/.bashrc", "-echo" );
                                [~, ~] = system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte library local installation setup <<<' >> ~/.bashrc", "-echo" );
                                [~, ~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                                [~, ~] = system( "chmod 777 ~/.bashrc && sh ~/.bashrc", "-echo" );
                            end
                            self.Err.msg    = "Whenever you intend to use ParaMonte in the future, before opening your MATLAB session, " + newline ...
                                            + "please execute the following command in your Bash shell to ensure all required paths " + newline ...
                                            + "are properly defined in your environment:" + newline + newline ...
                                            ... + "    " + string(strrep(setupFilePathCmd,'\','\\')) + newline + newline ...
                                            + "    " + setupFilePathCmd + newline + newline ...
                                            + "mission accomplished.";
                            self.Err.warn();
                        end

                    end

                    self.writeVerificationStatusFile("True");

                else

                    self.Err.msg = "Aborting the ParaMonte-for-MATLAB local build on your system.";
                    self.Err.warn();
                    return

                end

            end

            succeeded = true;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function warnForUnsupportedPlatform(self)
            self.Err.msg    = "The ParaMonte sampler kernel is currently exclusively available " ...
                            + "on AMD64 (64-bit) architecture for Windows/Linux/Darwin Operating Systems (OS). " ...
                            + "Your system appears to be of a different architecture or OS. As a result, " ...
                            + "the core sampler routines of ParaMonte will not be available on your system. However, " ...
                            + "the generic MATLAB interface of ParaMonte will is available on your system, which can " ...
                            + "be used for post-processing and visualization of the output files from " ...
                            + "already-performed ParaMonte simulations or other similar Monte Carlo simulations. " ...
                            + "There are ongoing efforts, right now as you read this message, to further increase the " ...
                            + "availability of ParaMonte library on a wider-variety of platforms and architectures. " ...
                            + "Stay tuned for updates by visiting " + newline + newline ...
                            + "    " + href(self.website.home.url) + newline + newline ...
                            + "That said," + newline + newline ...
                            + "if your platform is non-Windows and is compatible with GNU Compiler Collection (GCC)," + newline + newline ...
                            + "you can also build the required ParaMonte kernel's shared object files on your system " ...
                            + "by calling ParaMonte module's build() function from within your MATLAB environment.";
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;
            self.Err.prefix = self.names.paramonte;
            self.Err.warn();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function setupUnixPath(self)

            bashrcContents = getBashrcContents();
            dlibcmd = "export LD_LIBRARY_PATH=" + self.path.lib + ":$LD_LIBRARY_PATH";
            if ~contains(bashrcContents,dlibcmd)
                [~,~] = system( "chmod 777 ~/.bashrc", "-echo");
                [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte shared library setup >>>' >> ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte shared library setup <<<' >> ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                [~,~] = system( "chmod 777 ~/.bashrc && sh ~/.bashrc", "-echo" );
            end

            localInstallDir = self.getLocalInstallDir();
            if ~isempty(localInstallDir.root)

                pathcmd = [];
                dlibcmd = [];
                if ~isempty(localInstallDir.gnu.bin); pathcmd = "export PATH=" + localInstallDir.gnu.bin + ":$PATH"; end
                if ~isempty(localInstallDir.gnu.lib); dlibcmd = "export LD_LIBRARY_PATH=" + localInstallDir.gnu.lib + ":$LD_LIBRARY_PATH"; end
                pathcmdIsNotEmpty = ~isempty(pathcmd);
                dlibcmdIsNotEmpty = ~isempty(dlibcmd);
                pathcmdNotInBashrcContents = ~contains(bashrcContents,pathcmd);
                dlibcmdNotInBashrcContents = ~contains(bashrcContents,dlibcmd);
                if pathcmdIsNotEmpty || dlibcmdIsNotEmpty
                    if pathcmdNotInBashrcContents || dlibcmdNotInBashrcContents
                        [~,~] = system( "chmod 777 ~/.bashrc", "-echo");
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local GNU installation setup >>>' >> ~/.bashrc", "-echo" );
                        if pathcmdIsNotEmpty
                            if pathcmdNotInBashrcContents
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc", "-echo" );
                            end
                        end
                        if dlibcmdIsNotEmpty
                            if dlibcmdNotInBashrcContents
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc", "-echo" );
                            end
                        end
                    end
                    if pathcmdNotInBashrcContents || dlibcmdNotInBashrcContents
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte local GNU installation setup <<<' >> ~/.bashrc", "-echo" );
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                        [~,~] = system( "chmod 777 ~/.bashrc && sh ~/.bashrc", "-echo" );
                    end
                end

                pathcmd = [];
                dlibcmd = [];
                if ~isempty(localInstallDir.mpi.bin); pathcmd = "export PATH=" + localInstallDir.mpi.bin + ":$PATH"; end
                if ~isempty(localInstallDir.mpi.lib); dlibcmd = "export LD_LIBRARY_PATH=" + localInstallDir.mpi.lib + ":$LD_LIBRARY_PATH"; end
                pathcmdIsNotEmpty = ~isempty(pathcmd);
                dlibcmdIsNotEmpty = ~isempty(dlibcmd);
                pathcmdNotInBashrcContents = ~contains(bashrcContents,pathcmd);
                dlibcmdNotInBashrcContents = ~contains(bashrcContents,dlibcmd);
                if pathcmdIsNotEmpty || dlibcmdIsNotEmpty
                    if pathcmdNotInBashrcContents || dlibcmdNotInBashrcContents
                        [~,~] = system( "chmod 777 ~/.bashrc", "-echo");
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo");
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local MPI installation setup >>>' >> ~/.bashrc", "-echo" );
                        if pathcmdIsNotEmpty
                            if pathcmdNotInBashrcContents
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc", "-echo" );
                            end
                        end
                        if dlibcmdIsNotEmpty
                            if dlibcmdNotInBashrcContents
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc", "-echo" );
                                [~,~] = system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc", "-echo" );
                            end
                        end
                    end
                    if pathcmdNotInBashrcContents || dlibcmdNotInBashrcContents
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte local MPI installation setup <<<' >> ~/.bashrc", "-echo" );
                        [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                        [~,~] = system( "chmod 777 ~/.bashrc && sh ~/.bashrc", "-echo" );
                    end
                end
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function localInstallDir = getLocalInstallDir(self)

            localInstallDir = struct();
            localInstallDir.root = [];

            localInstallDir.mpi = struct();
            localInstallDir.mpi.root = [];
            localInstallDir.mpi.bin = [];
            localInstallDir.mpi.lib = [];

            localInstallDir.gnu = struct();
            localInstallDir.gnu.root = [];
            localInstallDir.gnu.bin = [];
            localInstallDir.gnu.lib = [];

            localInstallDir.caf = struct();
            localInstallDir.caf.root = [];
            localInstallDir.caf.bin = [];
            localInstallDir.caf.lib = [];

            pmGitRootDir = fullfile( self.path.lib , "paramonte-master" );

            if isfolder(pmGitRootDir)

                localInstallDir.root = pmGitRootDir;

                % mpi

                temp = fullfile( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "mpich", "3.2" );
                if isfolder(temp)
                    localInstallDir.mpi.root = temp;
                    temp = fullfile(localInstallDir.mpi.root, "bin");
                    if isfolder(temp); localInstallDir.mpi.bin = temp; end
                    temp = fullfile(localInstallDir.mpi.root, "lib");
                    if isfolder(temp); localInstallDir.mpi.lib = temp; end
                end

                % gnu

                temp = fullfile( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "gnu", "8.3.0" );
                if isfolder(temp)
                    localInstallDir.gnu.root = temp;
                    temp = fullfile(localInstallDir.gnu.root, "bin");
                    if isfolder(temp); localInstallDir.gnu.bin = temp; end
                    temp = fullfile(localInstallDir.gnu.root, "lib64");
                    if isfolder(temp); localInstallDir.gnu.lib = temp; end
                end

                % caf

                temp = fullfile( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "opencoarrays", "2.8.0" );
                if isfolder(temp)
                    localInstallDir.caf.root = temp;
                    temp = fullfile(localInstallDir.caf.root, "bin");
                    if isfolder(temp); localInstallDir.caf.bin = temp; end
                    temp = fullfile(localInstallDir.caf.root, "lib64");
                    if isfolder(temp); localInstallDir.caf.lib = temp; end
                end

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function mpi = findMPI(self)
            %
            % Return an struct containing the paths to different components of the MPI library.
            %
            mpi.path.broken = false;
            mpi.install.found = false;
            mpi.install.bin.found = false;
            mpi.install.bin.path = [];
            mpi.install.bin.mpiexec.found = false;
            mpi.install.bin.mpiexec.path = [];
            mpi.install.bin.mpivars.found = false;
            mpi.install.bin.mpivars.path = [];

            if self.platform.isWin32

                pathList = getenv("PATH");
                pathList = string(strsplit(pathList,';'));
                for thisPath = pathList

                    pathLower = lower(thisPath);
                    pathLower = string(strrep(pathLower,'\',''));
                    if contains(pathLower,"mpiintel64bin")

                        mpi.install.bin.found = true;
                        mpi.install.bin.path = thisPath;

                        mpiexecFilePath = string( fullfile(mpi.install.bin.path,"mpiexec.exe") );
                        if isfile(mpiexecFilePath)
                            mpi.install.bin.mpiexec.found = true;
                            mpi.install.bin.mpiexec.path = mpiexecFilePath;
                        end

                        mpivarsFilePath = string( fullfile(mpi.install.bin.path,"mpivars.bat") );
                        if isfile(mpivarsFilePath)

                            mpi.install.bin.mpivars.found = true;
                            mpi.install.bin.mpivars.path = mpivarsFilePath;

                            mpivarsCommand = """" + mpivarsFilePath + """";
                            self.Err.msg    = "Intel MPI library for 64-bit architecture detected at: " + newline ...
                                            + newline ...
                                            ... + "    " + string(strrep(mpi.install.bin.path,'\','\\')) + newline ...
                                            + "    """ + mpi.install.bin.path + """" + newline ...
                                            + newline ...
                                            + "To perform ParaMonte simulations in parallel on a single node, run the " + newline ...
                                            + "following two commands, in the form and order specified, on a MATLAB-aware " + newline ...
                                            + "mpiexec-aware command-line interface, such as the Intel Parallel Studio's command-prompt:" + newline ...
                                            + newline ...
                                            ... + "    " + string(strrep(mpivarsCommand,'\','\\')) + newline ...
                                            + "    " + mpivarsCommand + newline ...
                                            + newline ...
                                            + "    mpiexec -localonly -n 3 matlab -batch ""main_mpi""" + newline ...
                                            + newline ...
                                            + "where, " + newline ...
                                            + newline ...
                                            + "    0.   the first command defines the essential environment variables and, " + newline ...
                                            + "         the second command runs in the simulation in parallel, in which, " + newline ...
                                            + "    1.   you should replace the default number 3 with the number of " + newline ...
                                            + "         processors you wish to assign to your simulation task, " + newline ...
                                            + "    2.   main_mpi.m is the MATLAB file which serves as the entry point to " + newline ...
                                            + "         your simulation, where you call the ParaMonte sampler routines. " + newline ...
                                            + "         ATTN: Notice the MATLAB filename must appear without the file extension." + newline ...
                                            + "         ATTN: Notice the MATLAB filename must be enclosed with quotation marks." + newline ...
                                            + "    3.   -localonly indicates a parallel simulation on only a single node (this " + newline ...
                                            + "         flag will obviate the need for MPI library credentials registration). " + newline ...
                                            + "         For more information, visit: " + newline ...
                                            + newline ...
                                            + "             " + href(self.website.intel.mpi.windows.url) + newline ...
                                            + newline ...
                                            + "Note that the above two commands must be executed on a command-line that recognizes " + newline ...
                                            + "both MATLAB and mpiexec applications, such as the Intel Parallel Studio's command-prompt. " + newline ...
                                            + "For more information, in particular, on how to register to run Hydra services " + newline ...
                                            + "for multi-node simulations on Windows servers, visit: " + newline ...
                                            + newline ...
                                            + "    " + href(self.website.home.url);
                            self.Err.marginTop = 1;
                            self.Err.marginBot = 1;
                            self.Err.prefix = self.names.paramonte;
                            self.Err.note();

                            setupFilePath = fullfile(self.path.lib, "setup.bat");
                            fid = fopen(setupFilePath,"w");
                            fprintf(fid,"@echo off\n");
                            %fprintf(fid,"cd " + string(strrep(mpi.install.bin.path, '\', '\\')) + " && mpivars.bat quiet\n");
                            %fprintf(fid,"cd " + string(strrep(self.path.lib, '\', '\\')) + "\n");
                            fprintf(fid, "%s\n", "cd " + mpi.install.bin.path + " && mpivars.bat quiet");
                            fprintf(fid, "%s\n", "cd " + self.path.lib);
                            fprintf(fid,"@echo on\n");
                            fclose(fid);

                        end

                    end

                    mpi.install.found = mpi.install.bin.found && mpi.install.bin.mpiexec.found && mpi.install.bin.mpivars.found;
                    if mpi.install.found; break; end

                end

            elseif self.platform.isLinux

                pathList = getenv("PATH");
                pathList = string(strsplit(pathList,":"));
                for thisPath = pathList

                    pathLower = lower(thisPath);
                    pathLower = strrep(pathLower,"/","");
                    if contains(pathLower, "linuxmpiintel64")

                        mpi.install.bin.found = true;
                        mpi.install.bin.path = thisPath;

                        mpiexecFilePath = string( fullfile(mpi.install.bin.path,"mpiexec") );
                        if isfile(mpiexecFilePath)
                            mpi.install.bin.mpiexec.found = true;
                            mpi.install.bin.mpiexec.path = mpiexecFilePath;
                        end

                        mpivarsFilePath = string(fullfile(mpi.install.bin.path, "mpivars.sh"));

                        if isfile(mpivarsFilePath)

                            mpi.install.bin.mpivars.found = true;
                            mpi.install.bin.mpivars.path = mpivarsFilePath;

                            mpivarsCommand = """" + mpi.install.bin.mpivars.path + """";
                            self.Err.msg    = "Intel MPI library for 64-bit architecture detected at: " + newline + newline ...
                                            ... + "    " + string(strrep(mpi.install.bin.path,'\','\\')) + newline + newline ...
                                            + "    """ + mpi.install.bin.path + """" + newline + newline ...
                                            + "To perform ParaMonte simulations in parallel on a single node, run the " + newline ...
                                            + "following two commands, in the form and order specified, in a Bash shell, " + newline + newline ...
                                            ... + "    source " + string(strrep(mpivarsCommand,'\','\\')) + newline + newline ...
                                            + "    source " + mpivarsCommand + newline + newline ...
                                            + "    mpiexec -n 3 matlab -batch ""main_mpi""" + newline + newline ...
                                            + "where, " + newline + newline ...
                                            + "    0.   the first command defines the essential environment variables and, " + newline ...
                                            + "         the second command runs in the simulation in parallel, in which, " + newline ...
                                            + "    1.   you should replace the default number 3 with the number of " + newline ...
                                            + "         processors you wish to assign to your simulation task, " + newline ...
                                            + "    2.   main_mpi.m is the MATLAB file which serves as the entry point to " + newline ...
                                            + "         your simulation, where you call the ParaMonte sampler routines. " + newline ...
                                            + "         ATTN: Notice the MATLAB filename must appear without the file extension." + newline ...
                                            + "         ATTN: Notice the MATLAB filename must be enclosed with quotation marks." + newline + newline ...
                                            + "For more information on how to install and use and run parallel ParaMonte " + newline ...
                                            + "simulations on Linux systems, visit:" + newline + newline ...
                                            + "    " + href(self.website.home.url);
                            self.Err.marginTop = 1;
                            self.Err.marginBot = 1;
                            self.Err.prefix = self.names.paramonte;
                            self.Err.note();

                            try
                                setupFilePath = fullfile(self.path.lib, "setup.sh");
                                fid = fopen(setupFilePath,"w");
                                %fprintf(fid,"source " + string(strrep(mpivarsCommand, '\', '\\')));
                                fprintf(fid, "%s\n", "source " + mpivarsCommand);
                                fclose(fid);
                            catch
                                self.Err.msg    = "Failed to create the MPI setup file. " + newline ...
                                                + "It looks like the ParaMonte library directory is read-only. " + newline ...
                                                + "This can be potentially problematic. Skipping for now...";
                                self.Err.warn();
                            end

                        end

                    end

                    mpi.install.found = mpi.install.bin.found && mpi.install.bin.mpiexec.found && mpi.install.bin.mpivars.found;
                    if mpi.install.found; break; end

                end

            elseif self.platform.isMacOS

                gfortranPath = [];
                try
                    [~,gfortranVersion] = system("gfortran --version");
                    if contains(string(gfortranVersion), "GCC 9.")
                        [~,gfortranPath] = system("command -v gfortran");
                    end
                catch
                    warning("Failed to capture the GNU Compiler Collection version...");
                end

                mpi.install.bin.mpiexec.path = [];
                try
                    [~,mpiexecVersion] = system("mpiexec --version");
                    if contains(string(mpiexecVersion), "open-mpi")
                        [~,mpi.install.bin.mpiexec.path] = system("command -v mpiexec");
                    end
                catch
                    warning("Failed to capture the mpiexec version...");
                end

                if ~isempty(mpi.install.bin.mpiexec.path) && ~isempty(gfortranPath)
                    mpi.install.bin.found = true;
                    mpi.install.bin.mpiexec.found = true;
                    mpi.install.bin.mpivars.found = true; % dummy
                    [mpi.install.bin.path,~,~] = fileparts(mpi.install.bin.mpiexec.path); mpi.install.bin.path = string(mpi.install.bin.path);
                    self.Err.msg    = "MPI runtime libraries detected at: " + newline + newline ...
                                    ... + "    " + string(strrep(mpi.install.bin.path,'\','\\')) + newline + newline ...
                                    + "    " + mpi.install.bin.path + newline + newline ...
                                    + "To perform ParaMonte simulations in parallel on a single node, run the " + newline ...
                                    + "following command, in the form and order specified, in a Bash shell, " + newline + newline ...
                                    + "    mpiexec -n 3 matlab -batch ""main_mpi""" + newline + newline ...
                                    + "where, " + newline + newline ...
                                    + "    0.   the first command defines the essential environment variables and, " + newline ...
                                    + "         the second command runs in the simulation in parallel, in which, " + newline ...
                                    + "    1.   you should replace the default number 3 with the number of " + newline ...
                                    + "         processors you wish to assign to your simulation task, " + newline ...
                                    + "    2.   main_mpi.m is the MATLAB file which serves as the entry point to " + newline ...
                                    + "         your simulation, where you call the ParaMonte sampler routines. " + newline ...
                                    + "         ATTN: Notice the MATLAB filename must appear without the file extension." + newline ...
                                    + "         ATTN: Notice the MATLAB filename must be enclosed with quotation marks." + newline + newline ...
                                    + "For more information on how to install and use and run parallel ParaMonte " + newline ...
                                    + "simulations on Darwin operating systems, visit:" + newline + newline ...
                                    + "    " + href(self.website.home.url);
                    self.Err.marginTop = 1;
                    self.Err.marginBot = 1;
                    self.Err.prefix = self.names.paramonte;
                    self.Err.note();
                end

                mpi.install.found = mpi.install.bin.found && mpi.install.bin.mpiexec.found && mpi.install.bin.mpivars.found;

            else

                LocalInstallDir = getLocalInstallDir();
                if ~isempty(LocalInstallDir.mpi.bin) && ~isempty(LocalInstallDir.mpi.lib)

                    mpi.install.bin.found = true;
                    mpi.install.bin.path = LocalInstallDir.mpi.bin;

                    mpiexecFilePath = string( fullfile(mpi.install.bin.path,"mpiexec") );
                    if isfile(mpiexecFilePath)
                        mpi.install.bin.mpiexec.found = true;
                        mpi.install.bin.mpiexec.path = mpiexecFilePath;
                    end

                end

                mpi.install.bin.mpivars.found = mpi.install.bin.found && mpi.install.bin.mpiexec.found; % dummy
                mpi.install.found = mpi.install.bin.found && mpi.install.bin.mpiexec.found && mpi.install.bin.mpivars.found;

            end

            % one last try to find MPI library if not found

            if ~mpi.install.found

                mpi.path.broken = true;

                if self.platform.isLinux

                    defaultIntelLinuxMpiPath = self.getDefaultIntelLinuxMpiPath();
                    if defaultIntelLinuxMpiPath.mpiRootDirNotFound
                        return
                    else
                        mpi.install.found = true;
                        self.Err.msg    = "The PATH environmental variable of your Bash terminal does not point to " + newline ...
                                        + "any current installation of the Intel MPI runtime libraries on your system, " + newline ...
                                        + "however, ParaMonte was able to detect a hidden installation of Intel MPI runtime " + newline ...
                                        + "libraries on your system at, " + newline ...
                                        + newline ...
                                        + "    " + string(defaultIntelLinuxMpiPath.mpiDefaultRootDirList(end)) + newline ...
                                        + newline ...
                                        + "Follow the instructions below to ensure that your system's " + newline ...
                                        + "MPI runtime libraries will be properly detected in the future.";
                        self.Err.warn();
                        mpi.install.bin.path = self.setupIntelLinuxMpiPath(defaultIntelLinuxMpiPath);
                    end

                elseif self.platform.isWin32

                    self.Err.msg    = "Failed to detect the Intel MPI library for 64-bit architecture." + newline ...
                                    + "Now searching through the installed applications..." + newline ...
                                    + "This may take some time.";
                    self.Err.marginTop = 1;
                    self.Err.marginBot = 1;
                    self.Err.prefix = self.names.paramonte;
                    self.Err.warn();
                    [errorOccurred, installedApp] = system("wmic product get Name, Version");
                    if errorOccurred
                        return
                    elseif contains(installedApp,"Intel MPI")
                        mpi.install.found = true;
                        self.Err.msg = "Possible Intel MPI installation detected:";
                        self.Err.marginTop = 0;
                        self.Err.marginBot = 1;
                        self.Err.prefix = self.names.paramonte;
                        self.Err.note();
                        self.Err.marginBot = 0;
                        installedAppList = strsplit(installedApp,newline);
                        for i = 1:length(installedAppList)
                            app = installedAppList{i};
                            app = strrep(app,char(13),''); % remove cr
                            app = strrep(app,char(10),''); % remove newline
                            if contains(app,"Intel MPI")
                                self.Err.msg = "    " + app;
                                self.Err.note();
                            end
                        end
                        self.Err.marginTop = 1;
                        self.Err.marginBot = 1;
                    end

                end

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function dependencyList = getDependencyList(self)
            fileName = ".dependencies_";
            if self.platform.isWin32; fileName = fileName + "windows"; end
            if self.platform.isLinux; fileName = fileName + "linux"; end
            if self.platform.isMacOS; fileName = fileName + "macos"; end
            fullFilePath = string(fullfile(self.path.auxil,fileName));
            if isfile(fullFilePath)
                contents = fileread(fullFilePath);
                dependencyList = [string(strtrim(strsplit(contents,newline)))];
                dependencyList = [dependencyList(~contains(dependencyList,"!"))]; % remove comment lines
            else
                dependencyList = [];
                if self.platform.isWin32 || self.platform.isLinux
                    self.Err.msg    = "The Integrity of the ParaMonte library appears to have been compromised. " + newline ...
                                    + "The following file is missing in the ParaMonte library: " + newline + newline ...
                                    + "    " + fullFilePath + newline + newline ...
                                    + "To fix this issue, download a fresh clone of the ParaMonte MATLAB library from, " + newline + newline ...
                                    + "    " + href(self.website.github.release.url) ...
                                    ;
                    self.Err.abort();
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function prereqs = getPrereqs(self, downloadEnabled)
            %
            % optional DependencyList must be a string array (not a cell array).
            %
            prereqs = struct();
            prereqs.mpi = struct();
            prereqs.mpi.intel = struct();
            prereqs.list = self.getDependencyList();

            if nargin==1
                downloadEnabled = true;
            end

            if self.platform.isLinux
                intelMpiFilePrefix = "l_mpi-rt_";
                intelMpiFileSuffix = ".tgz";
            elseif self.platform.isWin32
                intelMpiFilePrefix = "w_mpi-rt_p_";
                intelMpiFileSuffix = ".exe";
            else
                return
            end

            self.Err.prefix = self.names.paramonte;
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;

            for dependency = prereqs.list

                fullFilePath = fullfile( self.path.lib, dependency );

                if downloadEnabled

                    % search all releases for the dependency

                    thisVersion = self.version.kernel.dump();
                    while getVecLen(thisVersion)
                        try
                            fullFilePath = websave(fullFilePath, self.website.github.release.url + "/download/" + thisVersion + "/" + dependency);
                            break
                        catch
                            thisVersion = getPreviousVersion(thisVersion);
                        end
                    end

                    if ~getVecLen(thisVersion)
                        self.Err.msg    = "Exhausted all releases of the ParaMonte library in search " + newline ...
                                        + "of the prerequisites, but could not find: " + dependency + newline ...
                                        + "Please report this issue at " + newline ...
                                        + newline ...
                                        + "    " + href(self.website.github.issues.url) + newline ...
                                        + newline ...
                                        + "In the meantime, visit, " + newline ...
                                        + newline ...
                                        + "    " + href(self.website.home.url) + newline ...
                                        + newline ...
                                        + "for instructions to manually install the MPI library on your " + newline ...
                                        + "system. Aborting the automatic MPI installation by ParaMonte...";
                        self.Err.warn();
                        return
                    end

                end

                if contains(dependency,intelMpiFilePrefix) && contains(dependency,intelMpiFileSuffix)
                    prereqs.mpi.intel.fullFileName = string( dependency );
                    prereqs.mpi.intel.fullFilePath = string( fullFilePath );
                    prereqs.mpi.intel.fileName = strsplit(prereqs.mpi.intel.fullFileName,intelMpiFileSuffix); prereqs.mpi.intel.fileName = prereqs.mpi.intel.fileName(1);
                    prereqs.mpi.intel.version = strsplit(prereqs.mpi.intel.fileName,intelMpiFilePrefix); prereqs.mpi.intel.version = prereqs.mpi.intel.version(2);
                    %self.prereqs.mpi.intel.version = string( dependency{1}(1:end-4) );
                    %mpiFileName = self.prereqs.mpi.intel.version;
                    %self.prereqs.mpi.intel.version = strsplit(self.prereqs.mpi.intel.version,"_");
                    %self.prereqs.mpi.intel.version = string(self.prereqs.mpi.intel.version(end));
                end

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function succeeded = installMPI(self)

            succeeded = false;

            self.Err.prefix = self.names.paramonte;
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;

            if self.platform.isWin32 || self.platform.isLinux

                self.Err.msg    = "Downloading the ParaMonte parallel library prerequisites..." + newline ...
                                + "Please make sure your firewall allows access to the Internet.";
                self.Err.note();

                self.prereqs = self.getPrereqs(); % download the prerequisites

                self.Err.msg    = "Installing the Intel MPI library for 64-bit architecture..." + newline ...
                                ...+ "file location: " + string(strrep(self.prereqs.mpi.intel.fullFilePath,'\','\\'));
                                + "file location: " + self.prereqs.mpi.intel.fullFilePath;
                self.Err.note();

                self.Err.msg    = "Please do not change the default installation location of the MPI library suggested by the installer." + newline ...
                                + "If you do change the default path, the onus will be on you to ensure the path to the " + newline ...
                                + "MPI runtime libraries exist in the environmental PATH variable of your session.";
                self.Err.warn();

                if self.platform.isWin32
                    [errorOccurred, output] = system("cmd /C " + self.prereqs.mpi.intel.fullFilePath + " &", "-echo");
                    if errorOccurred
                        self.Err.msg    = "Intel MPI library installation might have failed. Exit flag: " + string(errorOccurred) + newline ...
                                        + "command prompt output message: " + string(output) + newline ...
                                        + "The likely cause of this failure is that you did not open your current MATLAB session " ...
                                        + "as an administrator on your system. This automatic installation of the MPI libraries requires you to run MATLAB " ...
                                        + "with administrator privileges. To do so, " + newline + newline ...
                                        + "    1. Close your current MATLAB session." + newline ...
                                        + "    2. Find your MATLAB application again, but before opening it, right-click with your mouse on the application." + newline ...
                                        + "    3. Once you right-click, you should see an option like ""Run as administrator"". Click on this option." + newline ...
                                        + "    3. This wil open your MATLAB session with administrator privileges. Run this ParaMonte script again." + newline + newline ...
                                        + "As an alternative to the above solution, or if the aboe solution also fails, " ...
                                        + "you can open a Windows Command Prompt (cmd.exe) and type the following on the command prompt," + newline  + newline ...
                                        + "    """ + self.prereqs.mpi.intel.fullFilePath + """" + newline  + newline ...
                                        + "When doing this, make sure to include the quotation marks on the command line." + newline ...
                                        + "If you do not know how to open a Windows Command Prompt, click on the Window key on your keyboard and start typing," + newline  + newline ...
                                        + "    cmd" + newline  + newline ...
                                        + "or," + newline  + newline ...
                                        + "    command prompt" + newline  + newline ...
                                        + "Normally, the first application in the search reult is the Windows command prompt (cmd.exe). " ...
                                        + "Click on this application to open an instance of it and copy paste the above command on the terminal and press the Enter key." ...
                                        ;
                        self.Err.warn();
                        succeeded = false;
                        return
                    else
                        self.writeVerificationStatusFile("True");
                        self.Err.msg    = "The Intel MPI library installation appears to have succeeded. " + newline ...
                                        + "Now close your MATLAB environment and the command-line interface " + newline ...
                                        + "and reopen a new fresh command prompt.";
                        self.Err.note();
                    end
                end

                if self.platform.isLinux

                    %untar(self.prereqs.mpi.intel.fullFilePath,self.path.lib);
                    mpiExtractDir = fullfile(self.path.lib, self.prereqs.mpi.intel.fileName);
                    mpiExtractDirExists = isfolder(mpiExtractDir);
                    if mpiExtractDirExists
                        %[errorOccurred, output] = system("rm -rf " + mpiExtractDir, "-echo");
                        [status, ~, ~] = rmdir(mpiExtractDir, 's');
                        if status==1
                            mpiExtractDirExists = false;
                        end
                    end
                    if mpiExtractDirExists
                        errorOccurred = false;
                    else
                        disp("untarring " + self.prereqs.mpi.intel.fullFilePath + " at " + self.path.lib);
                        [errorOccurred, output] = system("cd " + self.path.lib + " && tar zxvf " + self.prereqs.mpi.intel.fullFileName, "-echo");
                    end
                    %try
                    %catch
                    if errorOccurred
                        self.Err.msg    = "Unzipping of Intel MPI runtime library tarball failed with the following message: " + newline ...
                                        + output + newline ...
                                        + "Make sure you have tar software installed on your system and try again.";
                        self.Err.abort();
                    end

                    self.Err.msg    = ..."If needed, use the following serial number when asked by the installer:" + newline + newline ...
                                    ...+ "    C44K-74BR9CFG" + newline + newline ...
                                    "When asked to choose the installation directory: if this is your personal computer, choose " + newline + newline ...
                                    + "    'install as root'" + newline + newline ...
                                    + "in the graphical user interface that appears in your session. " + newline + newline ...
                                    + "Otherwise, if you are using ParaMonte on a public server, " + newline ...
                                    + "for example, on a supercomputer, choose the third option:" + newline + newline ...
                                    + "   'install as current user'";
                    self.Err.note();

                    if self.isGUI
                        mpiInstallScriptName = "install_GUI.sh";
                    else
                        mpiInstallScriptName = "install.sh";
                    end
                    mpiInstallScriptPath = fullfile(mpiExtractDir, "install_GUI.sh");
                    if ~isfile(mpiInstallScriptPath)
                        self.Err. msg   = "Internal error occurred." + newline ...
                                        + "Failed to detect the Intel MPI installation Bash script." + newline ...
                                        + "Please report this issue at " + newline + newline ...
                                        + "    " + href(self.website.github.issues.url) + newline + newline ...
                                        + "Visit, " ...
                                        + "    " + href(self.website.home.url) + newline + newline ...
                                        + "for instructions to build the ParaMonte library on your system.";
                        self.Err.abort();
                    end

                    %[errorOccurred, ~] = system("cd " + mpiExtractDir + " && chmod +x ./" + mpiInstallScriptName, "-echo");
                    %if errorOccurred
                    %    self.Err.msg    = "The following action failed:" + newline + newline ...
                    %                    + "chmod +x " + string(strrep(mpiInstallScriptPath,'\','\\')) + newline + newline ...
                    %                    + "skipping...";
                    %    self.Err.warn();
                    %end

                    originalDir = cd(mpiExtractDir);
                    [errorOccurred, ~] = system("cd " + mpiExtractDir + " && ./" + mpiInstallScriptName + "", "-echo");
                    cd(originalDir);
                    if errorOccurred
                        self.Err.msg    = "Intel MPI runtime libraries installation for 64-bit " + newline ...
                                        + "architecture appears to have failed. If this failure " + newline ...
                                        + "is unexpected, please report this error at:" + newline ...
                                        + newline ...
                                        + "    " + href(self.website.github.issues.url) + newline ...
                                        + newline ...
                                        + "Visit, " + newline ...
                                        + newline ...
                                        + "    " + href(self.website.home.url) + newline ...
                                        + newline ...
                                        + "for instructions to build the ParaMonte library on your system.";
                        self.Err.warn();
                        succeeded = false;
                        return
                    end

                    self.Err.msg    = "Intel MPI runtime libraries installation for " + newline ...
                                    + "64-bit architecture appears to have succeeded." + newline ...
                                    + "Searching for the MPI runtime environment setup file...";
                    self.Err.note();

                    % setupFilePath = fullfile( self.path.lib, "setup.sh" );

                    defaultIntelLinuxMpiPath = self.getDefaultIntelLinuxMpiPath();
                    if defaultIntelLinuxMpiPath.mpiRootDirNotFound
                        self.Err.msg    = "Failed to detect the installation root path for Intel MPI runtime libraries " + newline ...
                                        + "for 64-bit architecture on your system. If you specified a different installation " + newline ...
                                        + "root path at the time of installation, please copy and paste it below. " + newline ...
                                        + "Note that the installation root path is part of the path that replaces: " + newline + newline ...
                                        + "    " + "opt" + newline + newline ...
                                        + "in the following path: " + newline + newline ...
                                        ... + "    " + string(strrep(fullfile("opt", defaultIntelLinuxMpiPath.mpiTrunkDir),'\','\\'));
                                        + "    " + string(fullfile("opt", defaultIntelLinuxMpiPath.mpiTrunkDir));
                        self.Err.warn();
                    end

%                    while true
%
%                        defaultIntelLinuxMpiPath = self.getDefaultIntelLinuxMpiPath();
%                        if defaultIntelLinuxMpiPath.mpiRootDirNotFound
%
%                            self.Err.msg    = "Failed to detect the installation root path for Intel MPI runtime libraries " + newline ...
%                                            + "for 64-bit architecture on your system. If you specified a different installation " + newline ...
%                                            + "root path at the time of installation, please copy and paste it below. " + newline ...
%                                            + "Note that the installation root path is part of the path that replaces: " + newline + newline ...
%                                            + "    " + "opt" + newline + newline ...
%                                            + "in the following path: " + newline + newline ...
%                                            ... + "    " + string(strrep(fullfile("opt", defaultIntelLinuxMpiPath.mpiTrunkDir),'\','\\'));
%                                            + "    " + string(fullfile("opt", defaultIntelLinuxMpiPath.mpiTrunkDir));
%                            self.Err.warn();
%
%                            msg = newline ...
%                                + "    Please type the root path to the MPI installation below and press ENTER." + newline ...
%                                + "    If you don't know the root path, simply press ENTER to quit:" + newline ...
%                                + newline ...
%                                );
%                            defaultIntelLinuxMpiPath.installationRootDirList = string(strtrim(input(msg,"s")));
%                            if getVecLen(defaultIntelLinuxMpiPath.installationRootDirList)
%                                if ntry > 1
%                                continue;
%                            else
%                                self.Err.msg = "Skipping the MPI runtime library environmental path setup...";
%                                self.Err.warn();
%                                break;
%                            end
%
%                        end
%                    end

                    if defaultIntelLinuxMpiPath.mpiRootDirNotFound

                        self.Err.msg    = "Failed to find the MPI runtime environment setup file on your system." + newline ...
                                        + "This is highly unusual. Normally, Intel MPI libraries must be installed" + newline ...
                                        + "in the following directory:" + newline + newline ...
                                        ... + "    " + string(strrep(defaultIntelLinuxMpiPath.mpiDefaultRootDirList(1),'\','\\')) + newline + newline ...
                                        + "    " + defaultIntelLinuxMpiPath.mpiDefaultRootDirList(1) + newline + newline ...
                                        + "or," + newline + newline ...
                                        ... + "    " + string(strrep(defaultIntelLinuxMpiPath.mpiDefaultRootDirList(2),'\','\\')) + newline + newline ...
                                        + "    " + defaultIntelLinuxMpiPath.mpiDefaultRootDirList(2) + newline + newline ...
                                        + "If you cannot manually find the Intel MPI installation directory," + newline ...
                                        + "it is likely that the installation might have somehow failed." + newline ...
                                        + "If you do find the installation directory, try to locate the" + newline ...
                                        + "'mpivars.sh' file which is normally installed in the following path:" + newline + newline ...
                                        ... + "    " + string(strrep(defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(1),'\','\\')) + newline + newline ...
                                        + "    " + defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(1) + newline + newline ...
                                        + "or," + newline + newline ...
                                        ... + "    " + string(strrep(defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(2),'\','\\')) + newline + newline ...
                                        + "    " + defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(2) + newline + newline ...
                                        + "Before attempting to run any parallel ParaMonte simulation, " + newline ...
                                        + "make sure you source this file, like the following:" + newline + newline ...
                                        ... + "    source " + string(strrep(defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(1),'\','\\')) + newline + newline ...
                                        + "    source " + defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(1) + newline + newline ...
                                        + "or," + newline + newline ...
                                        ... + "    source " + string(strrep(defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(2),'\','\\')) + newline + newline ...
                                        + "    source " + defaultIntelLinuxMpiPath.mpivarsDefaultFilePathList(2) + newline + newline ...
                                        + "where you will have to replace the path in the above with the " + newline ...
                                        + "correct path that you find on your system.";
                        self.Err.warn();

                    else

                        mpiBinDir = self.setupIntelLinuxMpiPath(defaultIntelLinuxMpiPath); % returns mpiBinDir which is not used

                    end

                end

                succeeded = true;

            elseif self.platform.isMacOS

                self.Err.msg    = "To use the ParaMonte kernel routines in parallel on macOS, " + newline ...
                                + "the Open-MPI library will have to be installed on your system. " + newline ...
                                + "Building the ParaMonte library prerequisites on your system...";
                self.Err.warn();
                succeeded = self.buildParaMontePrereqsForMac();

            else

                self.Err.msg    = "To use ParaMonte in parallel on this unknown Operating System, " + newline ...
                                + "ParaMonte needs to be built from scratch on your system. " + newline ...
                                + "Building ParaMonte library prerequisites on your system...";
                self.Err.warn();
                succeeded = self.build();

            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function mpiBinDir = setupIntelLinuxMpiPath(self,defaultIntelLinuxMpiPath)
            mpiBinDir = string(fullfile(defaultIntelLinuxMpiPath.mpiDefaultRootDirList(end), "bin"));
            mpiLibDir = string(fullfile(defaultIntelLinuxMpiPath.mpiDefaultRootDirList(end), "lib"));
            mpivarsFilePath = string(fullfile(mpiBinDir, "mpivars.sh"));
            if isfile(mpivarsFilePath)
                setupFilePath = fullfile( self.path.lib, "setup.sh" );
                try
                    fid = fopen(setupFilePath,"w");
                    %fprintf(fid,string(strrep(mpiBinDir,'\','\\'))+"\n");
                    %fprintf(fid,string(strrep(mpiLibDir,'\','\\'))+"\n");
                    %fprintf(fid,"source "+string(strrep(mpivarsFilePath,'\','\\')));
                    fprintf(fid, "%s\n", mpiBinDir);
                    fprintf(fid, "%s\n", mpiLibDir);
                    fprintf(fid, "%s\n", "source "+ mpivarsFilePath);
                    fclose(fid);
                catch
                    self.Err.msg    = "Failed to create the MPI setup file. " + newline ...
                                    + "It looks like the ParaMonte library directory is read-only. " + newline ...
                                    + "This can be potentially problematic. Skipping for now...";
                    self.Err.warn();
                end

                self.Err.msg    = "To ensure all MPI routine environmental variables " + newline ...
                                + "are properly load, source the following Bash script " + newline ...
                                + "in your Bash environment before calling mpiexec, like:" + newline + newline ...
                                ... + "    source " + string(strrep(mpivarsFilePath,'\','\\')) + newline + newline ...
                                + "    source " + mpivarsFilePath + newline + newline ...
                                + "Alternatively, ParaMonte can also automatically add" + newline ...
                                + "the required script to your '.bashrc' file, so that" + newline ...
                                + "all required MPI environmental variables are loaded" + newline ...
                                + "automatically before any ParaMonte usage from any" + newline ...
                                + "Bash command line on your system.";
                self.Err.note();

                bashrcContents = getBashrcContents();
                mpivarsFileCommand = "source " + mpivarsFilePath;
                if ~contains(bashrcContents,mpivarsFileCommand)

                    [~,~] = system( "chmod 777 ~/.bashrc", "-echo");
                    [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                    [~,~] = system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte MPI runtime library initialization >>>' >> ~/.bashrc", "-echo" );
                    [~,~] = system( "chmod 777 ~/.bashrc && echo '" + string(strrep(mpivarsFileCommand,'\','\\')) + "' >>  ~/.bashrc", "-echo" );
                    [~,~] = system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte MPI runtime library initialization <<<' >> ~/.bashrc", "-echo" );
                    [~,~] = system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc", "-echo" );
                    [~,~] = system( "chmod 777 ~/.bashrc && sh ~/.bashrc", "-echo" );

                end

                self.Err.msg    = "If you intend to run parallel simulations right now, you need to, " + newline + newline ...
                                + "    1.  quit your current MATLAB session, " + newline ...
                                + "    2.  close your current shell environment (and all Bash terminals, if possible)," + newline ...
                                + "    3.  open a new fresh Bash terminal (this is essential), " + newline ...
                                + "    4.  navigate to the root directory of your parallel ParaMonte simulation code, " + newline ...
                                + "    5.  invoke the MPI launcher to run your code in parallel." + newline + newline ...
                                + "For detailed instructions on running parallel ParaMonte simulations in MATLAB, " + newline ...
                                + "in particular, step 5 in the above, visit the ParaMonte documentation website, " + newline + newline ...
                                + "    " + href(self.website.home.url) ...
                                ;
                self.Err.note();

            else

                self.Err.msg    = "ParaMonte was able to detect an MPI library path on your system, however," + newline ...
                                + "the MPI installation appears to be corrupted. The required mpivars.sh " + newline ...
                                + "does not exist:" + newline + newline ...
                                ... + string(strrep(mpivarsFilePath,'\','\\'));
                                + mpivarsFilePath;
                self.Err.abort();

            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function mpiPath = getDefaultIntelLinuxMpiPath(self)
            mpiPath = struct();
            mpiPath.mpiDefaultRootDirList = [];
            mpiPath.mpiRootDirNotFound = true;
            mpiPath.mpivarsDefaultFilePathList = [];
            mpiPath.installationRootDirList = [ "/opt", self.path.home ];
            mpiPath.mpiTrunkDir = fullfile("intel", "compilers_and_libraries_" + self.prereqs.mpi.intel.version, "linux", "mpi", "intel64");
            for installationRootDir = mpiPath.installationRootDirList
                mpiPath.mpiDefaultRootDirList = [ mpiPath.mpiDefaultRootDirList, string(fullfile(installationRootDir, mpiPath.mpiTrunkDir)) ];
                mpiPath.mpivarsDefaultFilePathList = [ mpiPath.mpivarsDefaultFilePathList , fullfile(mpiPath.mpiDefaultRootDirList(end),"bin","mpivars.sh") ];
                if isfolder(mpiPath.mpiDefaultRootDirList(end))
                    mpiPath.mpiRootDirNotFound = false;
                    break
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function succeeded = buildParaMontePrereqsForMac(self)

            succeeded = false;

            self.Err.prefix = self.names.paramonte;
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;

            self.Err.msg = "Checking if Homebrew exists on your system...";
            self.Err.note();

            [~,brewPath] = system("command -v brew", "-echo");
            if ~isfile(brewPath)

                self.Err.msg = "Failed to detect Homebrew on your system. Installing Homebrew...";
                self.Err.note();

                %%%% install xcode-select

                [errorOccurred, ~] = system('xcode-select --install', "-echo");
                if errorOccurred
                    [errorOccurred, cmdout] = system("xcode-select -p");
                    if errorOccurred || ~isfolder( strrep( strrep(string(cmdout),char(10),"") , char(13) , "" ) )
                        self.Err.msg = self.getMacosInstallHelpMsg("xcode-select");
                        self.Err.warn();
                        return;
                    end
                end

                %%%% install Homebrew

                [errorOccurred1, ~] = system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"', "-echo");
                [errorOccurred2, ~] = system('brew --version', "-echo");
                if errorOccurred1 || errorOccurred2
                    self.Err.msg = self.getMacosInstallHelpMsg("Homebrew");
                    self.Err.warn();
                    return;
                end

            end

            %%%% cmake

            cmakeInstallationNeeded = false;
            [~,cmakePath] = system("command -v cmake", "-echo"); cmakePath = string(cmakePath);
            if ~isfile(cmakePath)

                cmakeInstallationNeeded = true;
                self.Err.msg = "cmake installation is missing on your system.";
                self.Err.note();

            else

                %self.Err.msg = "cmake installation detected at: " + string(strrep(cmakePath,'\','\\')) + "\n" + "Checking cmake version...";
                self.Err.msg = "cmake installation detected at: " + cmakePath + "\n" + "Checking cmake version...";
                self.Err.marginTop  = 0;
                self.Err.marginBot  = 0;
                self.Err.note();

                [errorOccurred,cmakeVersion] = system("cmake --version", "-echo");
                if errorOccurred

                    cmakeInstallationNeeded = true;
                    self.Err.msg = "Failed to detect the current cmake installation version. skipping...";
                    self.Err.note();

                else

                    cmakeVersion = strsplit(cmakeVersion," "); cmakeVersion = strsplit(cmakeVersion{3},"-"); cmakeVersion = string(cmakeVersion{1});
                    cmakeVersionList = string(strsplit(cmakeVersion,"."));

                    self.Err.msg = "current cmake version: " + cmakeVersion;
                    self.Err.note();

                    if str2double(cmakeVersionList{1})>=3 && str2double(cmakeVersionList{2})>=14
                        self.Err.msg = "cmake version is ParaMonte-compatible!";
                        self.Err.note();
                    else
                        cmakeInstallationNeeded = true;
                        self.Err.msg = "cmake version is NOT ParaMonte-compatible.";
                        self.Err.note();
                    end

                end

            end

            if cmakeInstallationNeeded

                self.Err.msg = "Installing cmake...";
                self.Err.note();

                [errorOccurred1,~] = system("brew install cmake", "-echo");
                [errorOccurred2,~] = system("brew link --overwrite cmake", "-echo");
                if errorOccurred1 || errorOccurred2
                    self.Err.msg = "cmake installation or linking failed.";
                    self.Err.marginTop = 1;
                    self.Err.marginBot = 1;
                    self.Err.warn();
                    return
                end

                [~,cmakeVersion] = system("cmake --version", "-echo");
                cmakeVersion = strsplit(cmakeVersion," "); cmakeVersion = strsplit(cmakeVersion{3},"-"); cmakeVersion = string(cmakeVersion{1});
                cmakeVersionList = string(strsplit(cmakeVersion,"."));

                self.Err.marginTop = 1;
                self.Err.marginBot = 1;
                if str2double(cmakeVersionList{1})>=3 && str2double(cmakeVersionList{2})>=14
                    self.Err.msg = "cmake installation succeeded.";
                    self.Err.note();
                else
                    self.Err.msg = self.getMacosInstallHelpMsg("cmake");
                    self.Err.warn();
                    return;
                end

            end

            %%%% gnu

            self.Err.msg = "Installing GNU Compiler Collection...";
            self.Err.note();

            [errorOccurred1,~] = system("brew install gcc@10", "-echo");
            [errorOccurred2,~] = system("brew link gcc@10", "-echo");

            if errorOccurred1 || errorOccurred2
                self.Err.msg = self.getMacosInstallHelpMsg("GNU Compiler Collection");
                self.Err.warn();
                return;
            end

            %%%% open-mpi

            self.Err.msg = "Installing Open-MPI...";
            self.Err.note();

            [errorOccurred1,~] = system("brew install open-mpi", "-echo");
            [errorOccurred2,~] = system("brew link open-mpi", "-echo");

            if errorOccurred1 || errorOccurred2
                self.Err.msg = self.getMacosInstallHelpMsg("Open-MPI");
                self.Err.warn();
                return;
            end

            succeeded = true;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function msg = getMacosInstallHelpMsg(self,app)
            if nargin<2; app = ""; end
            msg = "Failed to install and link the application " + app + " on your " + newline ...
                + "system. The " + app + " application is required to install and " + newline ...
                + "build the ParaMonte components and prerequisites. " + newline ...
                + newline ...
                + "If this is the first installation try, then, a second installation " + newline ...
                + "attempt might resolve the error. To do so, type the following " + newline ...
                + "commands on your MATLAB command line, " + newline ...
                + newline ...
                + "    pm = paramonte();" + newline ...
                + "    pm.verify();" + newline ...
                + newline ...
                + "Otherwise, you can install the application " + app + " manually " + newline ...
                + "on your system. The " + app + " installation is only a single " + newline ...
                + "command and takes only a few seconds to install. " + newline ...
                + "You can get the installation command from this page: " + newline ...
                + newline ...
                + "    " + self.website.home.install.macos.prereqs.cmd.url + newline ...
                + newline ...
                + "Once you have manually installed the missing component, retry " + newline ...
                + "the following commands on your MATLAB command line, " + newline ...
                + newline ...
                + "    pm = paramonte();" + newline ...
                + "    pm.verify();" + newline ...
                + newline ...
                + "skipping the installation for now..." ...
                ;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function writeVerificationStatusFile(self,verificationEnabled)
            fid = fopen(self.verificationStatusFilePath, "w");
            fprintf(fid,verificationEnabled);
            fclose(fid);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %function dispFinalMessage(self)
            %isYes = getUserResponse ( newline ...
            %                        + "    Shall we quit MATLAB now to perform the aforementioned tasks (y/n)? " ...
            %                        );
            %if isYes
            %    exit;
            %else
            %    self.Err.msg = "Continuing with ParaMonte at the risk of failing to run any parallel simulations...";
            %    self.Err.warn();
            %end
        %end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function displayParaMonteBanner(self)
            bannerFilePath = fullfile( self.path.auxil, ".ParaMonteBanner");
            versionLen = self.version.interface.dump();
            versionLen = length(versionLen{1});
            offset = fix( (versionLen - 4) / 2 );
            text = fileread(bannerFilePath);
            text = strrep(text, string(repmat(' ',1,offset))+"Version 0.0.0", "Version "+self.version.interface.dump);
            text = strrep(text, string(char(13)),"");
            disp(newline+text+newline);
            %fprintf(1,"\n");
            %lineList = string(strsplit(text,newline));
            %if self.isGUI % (self.platform.isWin32 || self.platform.isMacOS) &&
            %    for lineElement = lineList
            %        if contains(lineElement,"Version")
            %            line = strrep(lineElement, string(repmat(' ',1,offset))+"Version 0.0.0", "Version "+self.version.interface.dump);
            %        else
            %            line = lineElement;
            %        end
            %        fprintf(1,line);
            %    end
            %else
            %    newText = strrep(text, string(repmat(' ',1,offset))+"Version 0.0.0", "Version "+self.version.interface.dump);
            %    fprintf(1,newText);
            %end
            %fprintf(1,"\n\n");
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function result = helpme()
        %   Returns the help documentation for the object as a string
        %
        %   Parameters
        %   ----------
        %       None
        %
        %   Returns
        %   -------
        %       None
            result = help(paramonte);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Static)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Static, Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function pathSetupMsg = getPathSetupMsg()
            pathSetupMsg    = "    1.  quit your current MATLAB session, " + newline ...
                            + "    2.  open a new Bash terminal, or instead type the following on your current Bash command line, " + newline ...
                            + newline ...
                            + "            source ~/.bashrc" + newline ...
                            + newline ...
                            + "    3.  then, " + newline ...
                            + newline ...
                            + "            - for serial applications: call matlab from the same Bash command line as in the above, " + newline ...
                            + "            - for parallel applications: call your script form the Bash command line.";
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % methods (Static)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
