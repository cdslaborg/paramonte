%>  \brief
%>  This is the base class for the ParaMonte sampler routines.
%>
%>  \note
%>  See below for information on the attributes (properties).
%>
%>  \note
%>  See below for information on the methods.
%>
%>  \final{Sampler}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:43 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Sampler < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  \param[in] input   :   Optional scalar MATLAB string representing the path to an external
        %>                          input file that contains the simulation specification namelist(s).
        %>
        %>  \warning
        %>  USE THIS OPTIONAL ARGUMENT ONLY IF YOU KNOW WHAT YOU ARE DOING.
        %>  Specifying an input simulation file will force the ParaMonte sampler
        %>  to ignore all other simulation specifications set by the user via
        %>  the ``spec`` component of the sampler instance.
        %>
        input = "";
        %>
        %>   \param checked :   Optional scalar MATLAB logical.
        %>                      If set to ``true``, then the checked (verified) versions of the
        %>                      ParaMonte MATLAB shared libraries will be used for simulations.
        %>                      This means that most function calls and computations will be
        %>                      checked at runtime to ensure correctness, at the cost of
        %>                      noticeably slower runtime simulation speeds.<br>
        %>                      Set this option to ``true`` only for diagnostics, testing, and development.
        %>                      If set to empty, the sampler will first search for the non-checked version
        %>                      and if none found, will search for a checked version of the library.
        %>
        checked = [];
        %>
        %>  \param  mpiname :   Optional scalar MATLAB string representing the (vendor) name of the
        %>                      MPI runtime library. The following three library names are supported:<br>
        %>                          -   "impi"      :   The Intel MPI runtime library (Windows/Linux).
        %>                          -   "mpich"     :   The MPICH MPI runtime library (Linux/macOS).
        %>                          -   "openmpi"   :   The OpenMPI runtime library (Linux/macOS).<br>
        %>                      By default, you should set this argument to ``pm.lib.mpi.choice()``
        %>                      to enable the right preferred choice of MPI-enabled ParaMonte libraries.
        %>                      If for any reason, you want to use a non-preferred non-default MPI library,
        %>                      you can set the MPI library name explicitly as in the above list.<br>
        %>                      Beware that you must launch MATLAB with an ``mpiexec`` binary that
        %>                      comes from the same MPI library vendor specified for ``mpiname``.
        %>                      Note that the MPI-parallel ParaMonte shared libraries
        %>                      corresponding to the specified value for ``mpiname``
        %>                      may not exist in the ParaMonte package.
        %>                      For example,<br>
        %>                          -   On Windows, only ParaMonte library
        %>                              builds with Intel MPI are supported.<br>
        %>                          -   On Linux, the ParaMonte library builds with
        %>                              all MPI libraries listed above are supported,
        %>                              but there is no guarantee of their availability in
        %>                              the package. Only the ParaMonte library builds with
        %>                              Intel MPI (``"impi"``) are guaranteed to exist,
        %>                              unless user builds the ParaMonte MATLAB package
        %>                              from GitHub source for the desired MPI library.<br>
        %>                          -   On Darwin (macOS), only the ParaMonte shared library
        %>                              builds with MPICH and OpenMPI libraries are supported,
        %>                              but there is no guarantee of their availability in the
        %>                              package. Only the ParaMonte library builds with OpenMPI
        %>                              (``"openmpi"``) are guaranteed to exist, unless user
        %>                              builds the ParaMonte MATLAB package from GitHub
        %>                              source for the desired specified MPI library.<br>
        %>                      \warning                        
        %>                      USE THIS OPTIONAL ARGUMENT ONLY IF YOU KNOW ITS RAMIFICATIONS.
        %>                      This option should rarely be needed in any normal simulation..<br>
        %>                      The default value for ``mpiname`` is an empty string,
        %>                          implying no use of MPI-enabled ParaMonte library routines.
        %>
        mpiname = "";
        %>
        %>  \param  silent  :   Optional logical (Boolean) indicator which is ``false`` by default.
        %>                      If it is set to ``true``, it will silence all output postprocessing messages.
        %>
        %>  \note
        %>  Setting ``mpiname`` to a non-empty string
        %>  will automatically set ``silent`` to ``true``.
        %>
        silent = false;
        %>
        %>  \param  spec    :   A MATLAB structure containing all simulation specifications.
        %>                      All specifications are set to appropriate default values at runtime.
        %>                      Set the `spec` attributes to the desired values of your choice
        %>                      to override the default simulation specifications.
        %>
        spec = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        %>
        %>  \param  libspec :   Optional scalar MATLAB string containing a list of colon-(:)-separated
        %>                      configurations of the kernel shared library to be used for the sampling task.
        %>                      The most important of all is the compiler suite with which the library is built.
        %>
        %>  \warning
        %>  USE THIS OPTIONAL ARGUMENT ONLY IF YOU KNOW WHAT YOU ARE DOING.
        %>  The specified value for ``libspec`` is parsed internally to match
        %>  a corresponding library path name within the ParaMonte library.
        %>  The simulations will fail if no such directory can be found.
        %>
        libspec = "nocheck";
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        %>  \brief
        %>  Return a scalar object of class ``pm.sampling.Sampler``.
        %>
        %>  \details
        %>  This is the constructor of the class ``pm.sampling.Sampler``.
        %>  This class is not meant to be accessed by the end users.
        %>  It merely serves as the blueprint for the sampler subclasses
        %>  accessible to the end users.
        %>
        %>  \param[in]  method  :   The input scalar MATLAB string containing
        %>                          the name of the target ParaMonte sampler.
        %>
        %>  \return
        %>  `self`              :   The output scalar object of class ``pm.sampling.Sampler``.
        %>
        %>  \interface{Sampler}
        %>  \code{.m}
        %>
        %>      sampler = pm.sampling.Sampler();
        %>      sampler = pm.sampling.Sampler([]);
        %>      sampler = pm.sampling.Sampler(method);
        %>
        %>  \endcode
        %>  \final{Sampler}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 12:52 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = Sampler(method)
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