classdef Sampler < pm.matlab.Handle
    %
    %   This is the base class for the ParaMonte sampler routines.
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
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties(Access = public)
        %
        %   input
        %
        %       Optional scalar MATLAB string representing the path to an external
        %       input file that contains the simulation specification namelist(s).
        %
        %       \warning
        %
        %           USE THIS OPTIONAL ARGUMENT ONLY IF YOU KNOW WHAT YOU ARE DOING.
        %           Specifying an input simulation file will force the ParaMonte sampler
        %           to ignore all other simulation specifications set by the user via
        %           the ``spec`` component of the sampler instance.
        %
        input = "";
        %
        %   checked
        %
        %       Optional scalar MATLAB logical.
        %       If set to ``true``, then the checked (verified) versions of the
        %       ParaMonte MATLAB shared libraries will be used for simulations.
        %       This means that most function calls and computations will be
        %       checked at runtime to ensure correctness, at the cost of
        %       noticeably slower runtime simulation speeds.
        %
        %       Set this option to ``true`` only for diagnostics, testing, and development.
        %       If set to empty, the sampler will first search for the non-checked version
        %       and if none found, will search for a checked version of the library.
        %
        checked = [];
        %
        %   mpiname
        %
        %       Optional scalar MATLAB string representing the (vendor) name of the
        %       MPI runtime library. The following three library names are supported:
        %
        %           -   "impi"      :   The Intel MPI runtime library (Windows/Linux).
        %           -   "mpich"     :   The MPICH MPI runtime library (Linux/macOS).
        %           -   "openmpi"   :   The OpenMPI runtime library (Linux/macOS).
        %
        %       By default, you should set this argument to ``pm.lib.mpi.choice()``
        %       to enable the right preferred choice of MPI-enabled ParaMonte libraries.
        %       If for any reason, you want to use a non-preferred non-default MPI library,
        %       you can set the MPI library name explicitly as in the above list.
        %
        %       Beware that you must launch MATLAB with an ``mpiexec`` binary that
        %       comes from the same MPI library vendor specified for ``mpiname``.
        %       Note that the MPI-parallel ParaMonte shared libraries
        %       corresponding to the specified value for ``mpiname``
        %       may not exist in the ParaMonte package.
        %       For example,
        %
        %           -   On Windows, only ParaMonte library
        %               builds with Intel MPI are supported.
        %
        %           -   On Linux, the ParaMonte library builds with
        %               all MPI libraries listed above are supported,
        %               but there is no guarantee of their availability in
        %               the package. Only the ParaMonte library builds with
        %               Intel MPI (``"impi"``) are guaranteed to exist,
        %               unless user builds the ParaMonte MATLAB package
        %               from GitHub source for the desired MPI library.
        %
        %           -   On Darwin (macOS), only the ParaMonte shared library
        %               builds with MPICH and OpenMPI libraries are supported,
        %               but there is no guarantee of their availability in the
        %               package. Only the ParaMonte library builds with OpenMPI
        %               (``"openmpi"``) are guaranteed to exist, unless user
        %               builds the ParaMonte MATLAB package from GitHub
        %               source for the desired specified MPI library.
        %
        %       \warning
        %
        %           USE THIS OPTIONAL ARGUMENT ONLY IF YOU KNOW ITS RAMIFICATIONS.
        %           This option should rarely be needed in any normal simulation.
        %
        %       The default value for ``mpiname`` is an empty string,
        %       implying no use of MPI-enabled ParaMonte library routines.
        %
        mpiname = "";
        %
        %   silent
        %
        %       Optional logical (Boolean) indicator which is ``false`` by default.
        %       If it is set to ``true``, it will silence all output postprocessing messages.
        %
        %       \note
        %
        %           Setting ``mpiname`` to a non-empty string
        %           will automatically set ``silent`` to ``true``.
        %
        silent = false;
        %
        %   spec
        %
        %       A MATLAB structure containing all simulation specifications.
        %       All specifications are set to appropriate default values at runtime.
        %       Set the `spec` attributes to the desired values of your choice
        %       to override the default simulation specifications.
        %
        spec = [];
    end

    properties(Access = public, Hidden)
        mexname = "pm_sampling";
        libtype = "shared";
        partype = "serial";
        memtype = "heap";
        bldtype = ""; % build type
        clstype = ""; % compiler/linker suite
        name = "self";
        weblinks = [];
        matpath = []; % MATLAB paths.
        method = "";
        ndim = [];
        %njob = [];
        nml = [];
        %
        %   libspec
        %
        %       Optional scalar MATLAB string containing a list of colon-(:)-separated
        %       configurations of the kernel shared library to be used for the sampling task.
        %       The most important of all is the compiler suite with which the library is built.
        %
        %       \warning
        %
        %           USE THIS OPTIONAL ARGUMENT ONLY IF YOU KNOW WHAT YOU ARE DOING.
        %           The specified value for ``libspec`` is parsed internally to match
        %           a corresponding library path name within the ParaMonte library.
        %           The simulations will fail if no such directory can be found.
        %
        libspec = "nocheck";
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        function self = Sampler(method)
            %
            %   Return a scalar object of class ``pm.sampling.Sampler``.
            %
            %   This is the constructor of the class ``pm.sampling.Sampler``.
            %   This class is not meant to be accessed by the end users.
            %   It merely serves as the blueprint for the sampler subclasses
            %   accessible to the end users.
            %
            %   Parameters
            %   ----------
            %
            %       method
            %
            %           The input scalar MATLAB string containing
            %           the name of the target ParaMonte sampler.
            %
            %   Returns
            %   -------
            %
            %       self
            %
            %           The output scalar object of class ``pm.sampling.Sampler``.
            %
            %   Interface
            %   ---------
            %
            %       sampler = pm.sampling.Sampler();
            %       sampler = pm.sampling.Sampler([]);
            %       sampler = pm.sampling.Sampler(method);
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if  nargin < 1
                self.method = "sampler";
            else
                self.method = string(method);
            end
            self.weblinks = pm.lib.weblinks();
            if  self.method == "ParaDRAM"
                self.spec = pm.sampling.SpecDRAM(self.method);
            elseif self.method == "sampler"
                self.spec = pm.sampling.SpecBase(self.method);
            else
                error   ( newline ...
                        + "Unrecognized or unsupported ParaMonte sampler: " + self.method + newline ...
                        + newline ...
                        );
            end
            %filePath = mfilename("fullpath"); addpath(genpath(filePath),"-begin");
            self.name = inputname(1);
        end

        reportList = readReport(self, pattern)
        restartList = readRestart(self, pattern)
        chainList = readChain(self, pattern, sep)
        sampleList = readSample(self, pattern, sep)
        progressList = readProgress(self, pattern, sep)

    end

    methods(Hidden)
        finalize(self);
        run(self, getLogFunc, ndim)
        fileList = findfile(self, ftype, pattern);
    end

end