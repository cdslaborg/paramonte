%>  \brief
%>  This is the base class for the ParaMonte sampler routines.
%>
%>  \note
%>  See below for information on the attributes (properties).
%>
%>  \note
%>  See below for information on the methods.
%>
%>  \note
%>  Given the basic primitive nature of this class, it does not have a custom constructor.<br>
%>  Users should always use the subclasses of this base class.<br>
%>
%>  \final{Sampler}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:43 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef Sampler < pm.matlab.Handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public)
        %>
        %>  ``input``
        %>
        %>  The ``public`` scalar MATLAB string representing the path to an external
        %>  input file that contains the simulation specification namelist(s).<br>
        %>
        %>  \warning
        %>  <b>Use this optional argument only if you know what you are doing</b>.<br>
        %>  Specifying an input simulation file will force the ParaMonte sampler
        %>  to ignore all other simulation specifications set by the user via
        %>  the ``spec`` component of the sampler instance.<br>
        %>
        input = "";
        %>
        %>  ``checked``
        %>
        %>  The ``public`` scalar MATLAB logical.<br>
        %>  If set to ``true``, then the checked (verified) versions of the
        %>  ParaMonte MATLAB shared libraries will be used for simulations.<br>
        %>  This means that most function calls and computations will be
        %>  checked at runtime to ensure correctness, at the cost of
        %>  noticeably slower runtime simulation speeds.<br>
        %>  Set this option to ``true`` only for diagnostics, testing, and development.<br>
        %>  If set to empty, the sampler will first search for the non-checked version
        %>  and if none found, will search for a checked version of the library.<br>
        %>
        checked = [];
        %>
        %>  ``mpiname``
        %>
        %>  The ``public`` scalar MATLAB string representing the (vendor) name of the
        %>  MPI runtime library. The following three library names are supported:<br>
        %>  <ol>
        %>      <li>    ``"impi"``      :   The Intel MPI runtime library (Windows/Linux).<br>
        %>      <li>    ``"mpich"``     :   The MPICH MPI runtime library (Linux/macOS).<br>
        %>      <li>    ``"openmpi"``   :   The OpenMPI runtime library (Linux/macOS).<br>
        %>  </ol>
        %>  These names are normally one of the possible
        %>  values returned by [pm.lib.mpi.choices()](@ref choices).<br>
        %>  By default, the ParaMonte samplers should be capable of detecting
        %>  the right MPI runtime libraries and the ParaMonte library files.<br>
        %>  If the automated MPI detection fails or for any reason, or you want to use
        %>  a non-preferred non-default MPI library, you can set this argument explicitly
        %>  to one of the values returned by [pm.lib.mpi.choices()](@ref choices) to enable
        %>  your preferred choice of MPI-enabled ParaMonte libraries.<br>
        %>
        %>  Beware that you must launch MATLAB with an ``mpiexec`` binary that
        %>  comes with the same MPI library vendor string specified for ``mpiname``.<br>
        %>
        %>  The simulation task will mostly likely fail if the MPI-parallel ParaMonte
        %>  shared libraries corresponding to the specified value for ``mpiname``
        %>  do not exist in the ParaMonte package.<br>
        %>  For example,<br>
        %>  <ol>
        %>      <li>    On Windows, only ParaMonte library
        %>              builds with Intel MPI are generally supported.<br>
        %>      <li>    On Linux, the ParaMonte library builds with
        %>              all MPI libraries listed above are supported,
        %>              but there is no guarantee of their availability in
        %>              the package. Only the MPI-parallel ParaMonte library
        %>              names returned by [pm.lib.mpi.choices()](@ref choices)
        %>              are guaranteed to exist.<br>
        %>      <li>    On Darwin (macOS), the ParaMonte shared library builds
        %>              with MPICH and OpenMPI libraries are generally supported,
        %>              but there is no guarantee of the availability of all of them
        %>              in the ParaMonte package, unless the corresponding name is among
        %>              the names returned by [pm.lib.mpi.choices()](@ref choices).<br>
        %>  </ol>
        %>
        %>  \warning
        %>  <b>Use this optional argument only if you know its ramifications</b>.<br>
        %>  This option should rarely be needed in any normal MPI simulation.<br>
        %>  The default value for ``mpiname`` is an empty string, implying
        %>  no specific use of MPI-enabled ParaMonte library routines,
        %>  unless the ``mpiexec`` launcher with more than one process
        %>  is detected at runtime.<br>
        %>
        mpiname = "";
        %>
        %>  ``silent``
        %>
        %>  The ``public`` scalar MATLAB logical (Boolean) indicator which is ``false`` by default.<br>
        %>  If it is set to ``true``, it will silence all output postprocessing messages.<br>
        %>
        %>  \note
        %>  Setting ``mpiname`` to a non-empty string
        %>  will automatically set ``silent`` to ``true``.<br>
        %>
        silent = false;
        %>
        %>  ``spec``
        %>
        %>  The ``public`` MATLAB structure containing all simulation specifications.<br>
        %>  All specifications are set to appropriate default values at runtime.<br>
        %>  Set the `spec` attributes to the desired values of your choice
        %>  to override the default simulation specifications.<br>
        %>
        spec = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = public, Hidden)
        %>
        %>  ``mexname``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the name of the MEX file containing all ParaMonte sampler interfaces.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        mexname = "pm_sampling";
        %>
        %>  ``libtype``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the type of the ParaMonte Fortran library that is to be called.<br>
        %>  Although both ``static`` and ``shared`` may be possible, this value must always remain as ``shared``.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        libtype = "shared";
        %>
        %>  ``partype``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the parallelization type of the ParaMonte Fortran library that is to be called.<br>
        %>  Although the default value is ``serial``, the class methods will internally set this component to the right parallelization type at runtime.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        partype = "serial";
        %>
        %>  ``memtype``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the memory usage type of the ParaMonte Fortran library that is to be called.<br>
        %>  Although both ``heap`` and ``stack`` may be possible, this value must always remain as ``heap``.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        memtype = "heap";
        %>
        %>  ``bldtype``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the build type of the ParaMonte Fortran library that is to be called.<br>
        %>  The value of this component is determined at runtime depending on the availability of the ParaMonte Fortran library.<br>
        %>  The runtime value is among the values returned by the ParaMonte MATLAB library function [pm.lib.bldtypes()](@ref bldtypes).<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        bldtype = ""; % build type
        %>
        %>  ``clstype``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the type of the ParaMonte Fortran library parallelization that is to be called.<br>
        %>  Although the default value is ``serial``, the class methods will internally set this component to the right parallelization type at runtime.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        %>  \devnote
        %>  The attribute name ``clstype`` stands for *Compiler/Linker Suit* type.<br>
        %>
        clstype = ""; % compiler/linker suite
        %>
        %>  ``name``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the name of the
        %>  user-created object of superclass [pm.sampling.Sampler](@ref Sampler).<br>
        %>  This name is determined at runtime via the MATLAB intrinsic function ``inputname()``.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        name = "sampler";
        %>
        %>  ``weblinks``
        %>
        %>  A ``Hidden`` scalar MATLAB object of class [pm.lib.weblinks](@ref weblinks) containing the relevant ParaMonte library web links.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        weblinks = [];
        %>
        %>  ``matpath``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the contents of the MATLAB ``path`` variable at runtime.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        matpath = []; % MATLAB paths.
        %>
        %>  ``method``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing the name of the ParaMonte sampler subclass that the user has instantiated at runtime.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        method = "";
        %>
        %>  ``ndim``
        %>
        %>  A ``Hidden`` scalar MATLAB whole-number containing the user-specified value
        %>  for the number of dimensions of the mathematical density function to be sampled.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        ndim = [];
        %>
        %>  ``nml``
        %>
        %>  A ``Hidden`` scalar MATLAB string containing all the sampler specifications to be passed to the sampler as a Fortran namelist.<br>
        %>
        %>  \warning
        %>  This ``Hidden`` class component must not be accessed by the end users of the library.<br>
        %>
        nml = [];
        %>
        %>  ``libspec``
        %>
        %>  Optional scalar MATLAB string containing a list of colon-(``:``)-separated
        %>  configurations of the kernel shared library to be used for the sampling task.<br>
        %>  The most important of all is the compiler suite with which the library is built.<br>
        %>
        %>  \warning
        %>  <b>Use this optional argument only if you know what you are doing</b>.<br>
        %>  The specified value for ``libspec`` is parsed internally to match a
        %>  corresponding library path name within the ParaMonte library.<br>
        %>  The simulations will fail if no such directory can be found.<br>
        %>
        libspec = "nocheck";
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Return a scalar object of class [pm.sampling.Sampler``.<br>
        %>
        %>  \details
        %>  This is the constructor of the class [pm.sampling.Sampler](@ref Sampler).<br>
        %>  This class is not meant to be accessed by the end users.<br>
        %>  It merely serves as the blueprint for the subclasses of [pm.sampling.Sampler](@ref Sampler)
        %>  accessible to the end users.<br>
        %>
        %>  \param[in]  method  :   The input scalar MATLAB string containing
        %>                          the name of the target ParaMonte sampler.<br>
        %>
        %>  \return
        %>  ``self``            :   The output scalar object of class [pm.sampling.Sampler](@ref Sampler).<br>
        %>
        %>  \interface{Sampler}
        %>  \code{.m}
        %>
        %>      sampler = pm.sampling.Sampler();
        %>      sampler = pm.sampling.Sampler([]);
        %>      sampler = pm.sampling.Sampler(method);
        %>
        %>  \endcode
        %>
        %>  \final{Sampler}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 12:52 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
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
            elseif self.method == "ParaNest"
                self.spec = pm.sampling.SpecNest(self.method);
            elseif self.method == "sampler"
                self.spec = pm.sampling.SpecBase(self.method);
            else
                error   ( newline ...
                        + "Unrecognized or unsupported ParaMonte sampler: " + self.method + newline ...
                        + newline ...
                        );
            end
            %filePath = mfilename("fullpath"); addpath(genpath(filePath),"-begin");
            self.name = string(inputname(1));
        end

        reportList = readReport(self, pattern)
        restartList = readRestart(self, pattern)
        chainList = readChain(self, pattern, sep)
        sampleList = readSample(self, pattern, sep)
        progressList = readProgress(self, pattern, sep)
        run(self, getLogFunc, ndim)

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Hidden)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fileList = findfile(self, ftype, pattern);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ppm = getppm(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        finalize(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end