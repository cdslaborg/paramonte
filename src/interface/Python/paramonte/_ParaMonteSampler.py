####################################################################################################################################
####################################################################################################################################
####
####   MIT License
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   Permission is hereby granted, free of charge, to any person obtaining a 
####   copy of this software and associated documentation files (the "Software"), 
####   to deal in the Software without restriction, including without limitation 
####   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
####   and/or sell copies of the Software, and to permit persons to whom the 
####   Software is furnished to do so, subject to the following conditions:
####
####   The above copyright notice and this permission notice shall be 
####   included in all copies or substantial portions of the Software.
####
####   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
####   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
####   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
####   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
####   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
####   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
####   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
####
####   ACKNOWLEDGMENT
####
####   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
####   As per the ParaMonte library license agreement terms, if you use any parts of 
####   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
####   work (education/research/industry/development/...) by citing the ParaMonte 
####   library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import os
import sys
import typing as tp
import ctypes as ct

import _paramonte as pm
import _SpecBase as SpecBase
import _SpecMCMC as SpecMCMC
import _SpecDRAM as SpecDRAM
from _ReportFileContents import ReportFileContents
from _RestartFileContents import RestartFileContents
from _TabularFileContents import TabularFileContents

Struct = pm.Struct
newline = pm.newline

####################################################################################################################################
#### ParaMonteSampler class
####################################################################################################################################

class ParaMonteSampler:
    """

    This is the **ParaMonteSampler** base class for the ParaMonte
    sampler routines. This class is NOT meant to be directly accessed
    or called by the user of the ParaMonte library. However, its children,
    such the ParaDRAM sampler class will be directly accessible to the public.

        **Parameters**

            methodName

                A string representing the name of the ParaMonte sampler 
                that is to be instantiated.

        **Attributes**

            buildMode

                optional string argument with the default value "release".
                possible choices are:

                    "debug"

                        to be used for identifying sources of bug
                        and causes of code crash.

                    "release"

                        to be used in all other normal scenarios
                        for maximum runtime efficiency.

            mpiEnabled

                optional logical (boolean) indicator which is ``False`` by default.
                If it is set to ``True``, it will cause the ParaMonte simulation
                to run in parallel on the requested number of processors.
                See the class documentation guidelines in the above for
                information on how to run a simulation in parallel.

            reportEnabled

                optional logical (boolean) indicator which is ``True`` by default.
                If it is set to ``True``, it will cause extensive guidelines to be
                printed on the standard output as the simulation or post-processing
                continues with hints on the next possible steps that could be taken
                in the process. If you do not need such help and information set
                this variable to ``False`` to silence all output messages.

            inputFile

                optional string input representing the path to
                an external input namelist of simulation specifications.
                USE THIS OPTIONAL ARGUMENT WITH CAUTION AND
                ONLY IF YOU KNOW WHAT YOU ARE DOING.

                **WARNING**

                Specifying an input file will cause the sampler to ignore 
                all other simulation specifications set by the user via 
                sampler instance's `spec`-component attributes.

            spec

                A Python structure containing all simulation specifications.
                All simulation attributes are by default set to appropriate
                values at runtime. To override the default simulation
                specifications, set the ``spec`` attributes to some
                desired values of your choice.

                If you need help on any of the simulation specifications, try
                the supplied ``helpme()`` function in this component.

                If you wish to reset some specifications to the default values, 
                simply set them to ``None``.

        **Methods**

            See below for information on the methods.

        **Returns**

            Object of class ParaMonteSampler.

    """

    ################################################################################################################################
    #### ParaMonteSampler constructor
    ################################################################################################################################

    def __init__( self
                , methodName : str
                ):

        self._methodName = methodName
        self._libName = []
        self._ndim = []

        self._method = Struct()
        self._method.isParaDRAM = False
        self._method.isParaNest = False
        self._method.isParaTemp = False

        if self._methodName=="ParaDRAM":
            self._method.isParaDRAM = True
            self._objectName = "pmpd"
        else:
            pm.abort( msg   = "Internal error occurred. No sampling method other than ParaDRAM is currently " + newline
                            + "supported. Among its output files or simply, the path to the specific file " + newline
                            + "to be read. Please report this error at:" + newline
                            + newline
                            + "    " + pm.website.github.issues.url
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        self.buildMode = "release"
        self.mpiEnabled = False
        self.reportEnabled = not self.mpiEnabled
        self.inputFile = ""

        ############################################################################################################################
        #### ParaMonte specifications
        ############################################################################################################################

        self.spec = pm.utils.FrozenClass()

        # ParaMonte variables

        self.spec.sampleSize                               = None
        self.spec.randomSeed                               = None
        self.spec.description                              = None
        self.spec.outputFileName                           = None
        self.spec.outputDelimiter                          = None
        self.spec.chainFileFormat                          = None
        self.spec.variableNameList                         = None
        self.spec.restartFileFormat                        = None
        self.spec.outputColumnWidth                        = None
        self.spec.overwriteRequested                       = None
        self.spec.outputRealPrecision                      = None
        self.spec.silentModeRequested                      = None
        self.spec.domainLowerLimitVec                      = None
        self.spec.domainUpperLimitVec                      = None
        self.spec.parallelizationModel                     = None
        self.spec.progressReportPeriod                     = None
        self.spec.targetAcceptanceRate                     = None
        self.spec.mpiFinalizeRequested                     = None
        self.spec.maxNumDomainCheckToWarn                  = None
        self.spec.maxNumDomainCheckToStop                  = None

        if self._method.isParaDRAM:

            # ParaMCMC variables

            self.spec.chainSize                            = None
            self.spec.scaleFactor                          = None
            self.spec.startPointVec                        = None
            self.spec.proposalModel                        = None
            self.spec.proposalStartCovMat                  = None
            self.spec.proposalStartCorMat                  = None
            self.spec.proposalStartStdVec                  = None
            self.spec.sampleRefinementCount                = None
            self.spec.sampleRefinementMethod               = None
            self.spec.randomStartPointRequested            = None
            self.spec.randomStartPointDomainLowerLimitVec  = None
            self.spec.randomStartPointDomainUpperLimitVec  = None

            # ParaDRAM variables

            self.spec.adaptiveUpdateCount                  = None
            self.spec.adaptiveUpdatePeriod                 = None
            self.spec.greedyAdaptationCount                = None
            self.spec.delayedRejectionCount                = None
            self.spec.burninAdaptationMeasure              = None
            self.spec.delayedRejectionScaleFactorVec       = None

            self.spec.helpme = SpecDRAM.helpme

        self.spec._freeze()

    ################################################################################################################################
    #### _runSampler
    ################################################################################################################################

    def _runSampler ( self
                    , ndim          : int
                    , getLogFuncRaw : tp.Callable[[int,tp.List[float]], float]
                    , inputFile     : tp.Optional[str] = None
                    ) -> None:
        """

        Run ParaMonte sampler and return nothing. This method is identical to
        the ``runSampler()`` method, except that the input ``point`` parameter to
        the user-provided input objective function ``getLogFuncRaw(ndim,point)`` is
        a C-style raw pointer. This requires the user to guarantee that ``point`` will
        be always used with array bounds in their implementation of the objective function.
        The use of ``_runSampler()`` in place of ``runSampler()`` leads to a slight
        performance gain in the simulations.

            **Example serial usage**

            Copy and paste the following code enclosed between the
            two comment lines in your python/ipython/jupyter session
            (ensure the indentations of the pasted lines comply with Python rules):

            .. code-block:: python
                :linenos:

                ##################################
                import paramonte as pm
                import numpy as np
                def getLogFuncRaw(ndim,point):
                    # return the log of the standard multivariate
                    # Normal density function with ndim dimensions
                    return -0.5 * np.sum( np.double( point[0:ndim] )**2 )
                pmpd = pm.ParaDRAM()
                pmpd._runSampler( ndim = 4                      # length of point
                                , getLogFuncRaw = getLogFuncRaw # the objective function
                                )
                ##################################

            where,

                ndim

                    represents the number of dimensions of the domain of
                    the user's objective function ``getLogFuncRaw(ndim, point)``
                    and,

                getLogFuncRaw(ndim, point)

                    represents the user's objective function to be sampled,
                    where,

                        ndim

                            is a 32-bit integer, representing the number of
                            dimensions of the domain of the user-provided
                            objective function.

                        point

                            is a C-style array-pointer of length ``ndim``
                            and type float64. Note that the bounds of
                            ``point`` must be always specified wherever
                            it is used within the objective function.

                    On output, it must return the natural logarithm of
                    the objective function.

            **Parameters**

                All input parameters have the same meaning as the parameters
                of ``runSampler()``. The only difference is in the input
                parameters to the objective function ``getLogFuncRaw``.

            **Returns**

                None

        """

        #### verify ndim

        if not isinstance(ndim,int) or ndim<1:
            pm.abort( msg   = "The input argument ndim must be a positive integer," + newline
                            + "representing the number of dimensions of the domain of" + newline
                            + "the user's objective function ``getLogFunc()``." + newline
                            + "You have entered ndim = " + str(ndim)
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        #### verify getLogFuncRaw

        if not callable(getLogFuncRaw):
            pm.abort( msg   = "The input argument ``getLogFuncRaw`` must be a callable function." + newline
                            + "It represents the user's objective function to be sampled," + newline
                            + "which must take an input integer ndim representing the number of" + newline
                            + "dimensions of the domain of the objective function to be samples and," + newline
                            + "a second input argument of type numpy float64 array of length ndim." + newline
                            + "On return it must return the natural logarithm of the objective function."
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        #### verify mpiEnabled

        if not isinstance(self.mpiEnabled,bool):
            pm.abort( msg   = "The sampler attribute ``mpiEnabled`` must be of type bool." + newline
                            + "It is an optional logical (boolean) indicator which is False by default." + newline
                            + "If it is set to True, it will cause the ParaMonte simulation" + newline
                            + "to run in parallel on the requested number of processors." + newline
                            + "See the ParaMonte class information on how to run a simulation " + newline
                            + "in parallel. You have entered mpiEnabled = " + str(self.mpiEnabled)
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )
        if self.mpiEnabled: self.reportEnabled = False

        #### verify buildMode

        buildMode = "release" if self.buildMode is None else self.buildMode

        stype = None
        dummyList = None
        errorOccurred = True
        if isinstance(buildMode,str):
            errorOccurred = False
            if "-" in buildMode:
                dummyList = buildMode.split("-")
                buildMode = dummyList[0] # build type
                stype = dummyList[1] # compiler suite

        if not errorOccurred: errorOccurred = buildMode not in ["release", "testing", "debug"]
        if not errorOccurred and stype is not None: errorOccurred = not (stype=="gnu" or stype=="intel")
        if errorOccurred:
            if dummyList is not None: buildMode = "-".join(dummyList)
            pm.abort( msg   = "The object attribute ``buildMode`` must be of type ``str``." + newline
                            + "It is an optional string argument with default value \"release\"." + newline
                            + "possible choices are:" + newline
                            + newline
                            + "    \"debug\":" + newline
                            + newline
                            + "        to be used for identifying sources of bug" + newline
                            + "        and causes of code crash." + newline
                            + newline
                            + "    \"release\":" + newline
                            + newline
                            + "        to be used in all other normal scenarios" + newline
                            + "        for maximum runtime efficiency." + newline
                            + newline
                            + "You have entered buildMode = " + str(buildMode)
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        #### verify inputFile

        if inputFile is not None and not isinstance(inputFile,str):
            pm.abort( msg   = "The input argument ``inputFile`` must be of type str." + newline
                            + "It is an optional string input representing the path to" + newline
                            + "an external input namelist of simulation specifications." + newline
                            + "USE THIS OPTIONAL ARGUMENT WITH CAUTION AND" + newline
                            + "ONLY IF YOU KNOW WHAT YOU ARE DOING." + newline
                            + "Specifying this option will cause the ParaMonte sampler " + newline
                            + "to ignore all other simulation specifications set by " + newline
                            + "the user via sampler's instance ``spec`` attributes." + newline
                            + "You have entered inputFile = " + str(inputFile)
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        inputFileVec_pntr, inputFileLen = self._getInputFile(inputFile)

        if self.mpiEnabled:
            parallelism = "_mpi"
        else:
            parallelism = ""
            pm.note( msg   = "Running the " + self._methodName + " sampler in serial mode..." + newline
                            + "To run the " + self._methodName + " sampler in parallel mode visit:" + newline
                            + newline
                            + "    " + pm.website.home.url + newline
                            + newline
                            + "If you are using Jupyter notebook, check the Jupyter's " + newline
                            + "terminal window for realtime simulation progress and report."
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        #if len(sys.argv)>1:
        #    if sys.argv[1]=="p":
        #        pm.note( msg = Running sampler in parallel mode...
        #                , methodName = self._methodName
        #                )
        #        print("\nRunning sampler in parallel mode..." + newline)
        #        libName += "_mpi"
        #else:
        #    print("\nRunning ParaMonte sampler in serial mode..." + newline)
        #try:
        #    from mpi4py import MPI
        #    comm = MPI.COMM_WORLD
        #    libName += "_mpi"
        #    if comm.size==1:
        #        print("\nRunning ParaMonte sampler in serial mode..." + newline)
        #        if MPI.Is_initialized():
        #            print("Hello")
        #            MPI.Finalize()
        #    elif comm.rank==0:
        #        print("\nRunning ParaMonte sampler in parallel mode on {} processes..." + newline.format(comm.size))
        #        comm.barrier()
        #except ImportError:
        #    print("\nImportError occurred..." + newline)
        #    print("\nRunning ParaMonte sampler in serial mode..." + newline)

        sys.stdout.flush()

        #### setup env

        if pm.platform.isWin32:

            if "PATH" in os.environ:
                os.environ["PATH"] = pm.path.root + os.pathsep + os.environ["PATH"]
            else:
                os.environ["PATH"] = pm.path.root

            #mpiFound = False
            pathList = os.environ["PATH"].split(";")
            for path in pathList:
                pathLower = path.lower().replace("\\","")
                if ("mpiintel64bin" in pathLower):
                    #mpiFound = True
                    #bldMode = buildMode
                    #if bldMode=="testing": bldMode = "release"
                    mpiPath = os.path.join(path,"release")
                    os.environ["PATH"] = mpiPath + os.pathsep + os.environ["PATH"]
                    libfabricPath = os.path.join(os.path.dirname(path),"libfabric","bin")
                    os.environ["PATH"] = libfabricPath + os.pathsep + os.environ["PATH"]
                    break

        else:

            if "LD_LIBRARY_PATH" not in os.environ:
                os.environ["LD_LIBRARY_PATH"] = "."
                if self.reportEnabled and pm.platform.isLinux:
                    pm.warn( msg   = "LD_LIBRARY_PATH environmental variable is not defined in " + newline
                                    + "your Python session. Consider running the following command " + newline
                                    + "in your Bash shell before reopening Python and using ParaMonte: " + newline
                                    + newline
                                    + "    export LD_LIBRARY_PATH=."
                            , methodName = self._methodName
                            , marginTop = 1
                            , marginBot = 1
                            )
            libdir = "/usr/lib"
            if os.path.isdir(libdir):
                os.environ["LD_LIBRARY_PATH"]  = libdir + os.pathsep + os.environ["LD_LIBRARY_PATH"]
            libdir = "/usr/local/lib"
            if os.path.isdir(libdir):
                os.environ["LD_LIBRARY_PATH"]  = libdir + os.pathsep + os.environ["LD_LIBRARY_PATH"]
            libdir = "/usr/lib64"
            if os.path.isdir(libdir):
                os.environ["LD_LIBRARY_PATH"]  = libdir + os.pathsep + os.environ["LD_LIBRARY_PATH"]
            libdir = "/usr/local/lib64"
            if os.path.isdir(libdir):
                os.environ["LD_LIBRARY_PATH"]  = libdir + os.pathsep + os.environ["LD_LIBRARY_PATH"]
            os.environ["LD_LIBRARY_PATH"]  = pm.path.root + os.pathsep + os.environ["LD_LIBRARY_PATH"]

            from _pmreqs import getLocalInstallDir
            localInstallDir = getLocalInstallDir()

            if localInstallDir.gnu.root is not None:
                for object in os.scandir(localInstallDir.gnu.root):
                    if object.is_dir() and ("lib" in object.name):
                        os.environ["LD_LIBRARY_PATH"] = object.path + os.pathsep + os.environ["LD_LIBRARY_PATH"]

            if localInstallDir.mpi.root is not None:
                if localInstallDir.mpi.bin is not None: os.environ["PATH"] = localInstallDir.mpi.bin + os.pathsep + os.environ["PATH"]
                for object in os.scandir(localInstallDir.mpi.root):
                    if object.is_dir() and ("lib" in object.name):
                        os.environ["LD_LIBRARY_PATH"] = object.path + os.pathsep + os.environ["LD_LIBRARY_PATH"]
                if localInstallDir.mpi.lib is not None: os.environ["LD_LIBRARY_PATH"] = localInstallDir.mpi.lib + os.pathsep + os.environ["LD_LIBRARY_PATH"]

        # import ParaMonte dll define result (None) AND argument (pointer to a c function) type

        buildModeList = ["release","testing","debug"]
        buildModeList.pop(buildModeList.index(buildMode))
        buildModeList.insert(0,buildMode)
        pmcsList = ["intel","gnu"]
        if stype is not None:
            pmcsList.pop(pmcsList.index(stype))
            pmcsList.insert(0,stype)

        libNameSuffix = parallelism +   { "windows" : ".dll"
                                        , "cygwin"  : ".dll"
                                        , "darwin"  : ".dylib"
                                        , "linux"   : ".so"
                                        }.get(pm.platform.osname, ".so")

        libPath = None
        libFound = False
        libNamePrefix = "libparamonte_python_" + pm.platform.osname.lower() + "_" + pm.platform.arch + "_"
        from ctypes.util import find_library
        for buildMode in buildModeList:

            for pmcs in pmcsList:

                libName = libNamePrefix + pmcs + "_" + buildMode + "_dynamic_heap" + libNameSuffix
                libPath = find_library(libName)
                if libPath is None: libPath = os.path.join( pm.path.root, libName )

                libFound = os.path.isfile(libPath)
                if libFound: break

            if libFound: # check if lib file exists
                break
            #else:
            #    if self.reportEnabled:
            #        pm.warn( msg   = "ParaMonte dynamic library for the requested build mode " + buildMode + " not found." + newline
            #                        + "Searching for ParaMonte dynamic library in other build modes..."
            #                , methodName = self._methodName
            #                , marginTop = 1
            #                , marginBot = 1
            #                )
            #    #libName = libName.replace(buildMode,mode)
            #    #buildMode = mode

        if not libFound:
            if self.mpiEnabled:
                parallelMsg = ("This happens frequently with parallel simulations and the " + newline
                            + "most likely reason is that the user did NOT carefully follow " + newline
                            + "the ParaMonte instructions to successfully install and define " + newline
                            + "the variables of the MPI runtime library on their system. " + newline
                            + "To learn these about these instructions, type the following " + newline
                            + "in your Python session, " + newline
                            + newline
                            + "    import paramonte as pm" + newline
                            + "    pm.verify()" + newline
                            + newline
                            + "Then, carefully follow the instructions provided to define " + newline
                            + "the MPI runtime variables in your current Python session. " + newline
                            + "If the error still persists, please report this issue at: " + newline
                            )
            else:
                parallelMsg = "Please report this issue at:" + newline
            if pm.platform.isWin32:
                from _pmreqs import buildInstructionNoteWindows
                buildMsg = buildInstructionNoteWindows
            else:
                from _pmreqs import buildInstructionNoteUnix
                buildMsg = newline + newline + buildInstructionNoteUnix
            pm.abort( msg   = "Exhausted all possible ParaMonte dynamic library search" + newline
                            + "names but could not find any compatible library." + newline
                            #+ "Last search:" + newline
                            #+ newline
                            #+ "    " + libPath + newline
                            #+ newline
                            + "It appears your ParaMonte library is missing some files. " + newline
                            + parallelMsg
                            + newline
                            + "    " + pm.website.github.issues.url + newline
                            + newline
                            + "Visit," + newline
                            + newline
                            + "    " + pm.website.home.url + newline
                            + newline
                            + "for instructions on how to build the ParaMonte library" + newline
                            + "object files on your system."
                            + buildMsg
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        # define ctypes wrapper function, with the proper result and argument types
        _getLogFuncRaw_proc = ct.CFUNCTYPE ( ct.c_double                  # function result
                                           #, ct.POINTER(ct.c_int32)      # ndim
                                            , ct.c_int32                   # ndim
                                            , ct.POINTER(ct.c_double)     # point
                                            )
        getLogFuncRaw_pntr = _getLogFuncRaw_proc(getLogFuncRaw)

        try:

            pmdll = ct.CDLL(libPath)

        except Exception as e:

            import logging
            logger = logging.Logger("catch_all")
            logger.error(e, exc_info=True)

            from _pmreqs import buildInstructionNote
            pm.abort( msg  = "Failed to load the required ParaMonte shared library." + newline
                            + "This is either due to the incompatibility of the DLL with your" + newline
                            + "platform or due to missing some required dependent libraries." + newline
                            + "In either case, you can likely resolve this error by building." + newline
                            + "the required ParaMonte shared libraries on your system." + newline
                            + newline
                            + "Visit," + newline
                            + newline
                            + "    " + pm.website.home.url + newline
                            + newline
                            + "for instructions to build the ParaMonte library on your system." + newline
                            + newline
                            + "Please report this issue at:" + newline
                            + newline
                            + "    " + pm.website.github.issues.url + newline
                            + newline
                            + buildInstructionNote
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        pmdll.runParaDRAM.restype = ct.c_int32
        #pmdll.runParaDRAM.restype = None
        #pmdll.runParaDRAM.argtypes =    [ ct.POINTER(ct.c_int32)     # ndim
        pmdll.runParaDRAM.argtypes =    [ ct.c_int32                   # ndim
                                        , _getLogFuncRaw_proc           # procedure
                                        , ct.POINTER(ct.c_char)       # inputFile byte object
                                        , ct.c_int32                   # lenInpuFile
                                       #, ct.POINTER(ct.c_size_t)     # lenInpuFile
                                        , ]

        #def getLogFuncRawWrapper(ndim_pntr,point): return getLogFuncRaw(ndim[0],point)

        # construct procedure pointer
        #def getLogFuncRawWrapper(ndim,point): return getLogFuncRaw(np.array(point[0:ndim]))
        #getLogFuncRaw_pntr = _getLogFuncRaw_proc(getLogFuncRawWrapper)

        # construct ndim pointer
        #ndim_pntr = ct.byref(ct.c_int32(ndim))

        # call ParaMonte
        #pmdll.runParaDRAM   ( ndim_pntr
        #pmdll.runParaDRAM   ( ct.c_int32(ndim)


        if self._method.isParaDRAM:

            errFlag = pmdll.runParaDRAM ( ct.c_int32(ndim)
                                        , getLogFuncRaw_pntr
                                        , inputFileVec_pntr
                                        , ct.c_int32(inputFileLen)
                                       #, inputFileLen_pntr
                                        )

            if errFlag!=0:

                # first check for old existing files:

                existingFileList = []
                outputFileName = os.path.abspath(self.spec.outputFileName)
                if not os.path.isdir(outputFileName):
                    import glob
                    existingFileList = glob.glob(outputFileName+"*")

                if self._method.isParaDRAM and len(existingFileList)>1:
                    existingSimulationMsg = ( "It appears that an old " + self._methodName + " simulation with " + newline
                                            + "the same output file names exists in the specified path:" + newline
                                            + newline
                                            + newline.join(existingFileList)
                                            + newline + newline
                                            + "Keep in mind that old simulation files CANNOT BE OVERWRITTEN " + newline
                                            + "and attempting to do so will cause the new simulation to crash. " + newline
                                            + "If this is the case, specify a new filename prefix for the simulation " + newline
                                            + "output files or set it to `None` in the simulation's input specifications " + newline
                                            + "so that the sampler can automatically generate unique output " + newline
                                            + "filenames for your simulation."
                                            )
                else:
                    existingSimulationMsg = ""

                pm.abort( msg   = "The simulation failed. For more information, checkout the " + newline
                                + "simulation output error message on your your Bash / Python " + newline
                                + "terminal or command-prompt, also at the end of the output report " + newline
                                + "file, if it has been generated:" + newline
                                + newline
                                + "    " + outputFileName + "_report.txt" + newline
                                + newline
                                + existingSimulationMsg
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 1
                        )

        #def isLoaded(libPath):
        #    abslibPath =
        #    return os.system("lsof -p {} | grep {} > /dev/null".format( os.getpid(), os.path.abspath(libPath) )) == 0
        #def dlclose(libdll): libdll.dlclose(libdll._handle)

        if pm.platform.isWin32:
            handle = ct.windll.kernel32.LoadLibraryA(libPath)
            ct.windll.kernel32.FreeLibrary(handle)
        else:
           #while isLoaded(libPath):
           #    dlclose(pmdll._handle)
            try:
                ct.dlclose(pmdll._handle)
            except:
                if self.reportEnabled:
                    pm.warn ( msg   = "Failed to properly close the ParaMonte shared library file. " + newline
                                    + "This should not cause any major problems, unless you intend to " + newline
                                    + "run a new ParaMonte simulation, in which case, you may want to " + newline
                                    + "exit and re-enter your Python environment."
                            , methodName = self._methodName
                            , marginTop = 1
                            , marginBot = 1
                            )

        if self.reportEnabled:
            pm.note( msg    = "To read the generated output files, try:" + newline
                            + newline
                            + "    " + self._objectName + ".readReport()      # to read the summary report from the output report file." + newline
                            + "    " + self._objectName + ".readSample()      # to read the final i.i.d. sample from the output sample file." + newline
                            + "    " + self._objectName + ".readChain()       # to read the uniquely-accepted points from the output chain file." + newline
                            + "    " + self._objectName + ".readMarkovChain() # to read the Markov Chain. NOT recommended for very large chains." + newline
                            + "    " + self._objectName + ".readRestart()     # to read the contents of an ASCII-format output restart file." + newline
                            + "    " + self._objectName + ".readProgress()    # to read the contents of an output progress file." + newline
                            + newline
                            + "where you should replace `" + self._objectName + "` with your " + self._methodName + " sampler's object name." + newline
                            + "For more information and examples on the usage, visit:" + newline
                            + newline
                            + "    " + pm.website.home.url
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        return None

    ################################################################################################################################
    #### _getInputFile()
    ################################################################################################################################

    def _getInputFile(self, inputFile):

        if inputFile is None:

            ########################################################################################################################
            #### begin namelist generation from arguments
            ########################################################################################################################

            nameList = ""

            # setup outputFileName if it is None

            if self.spec.outputFileName is None:
                self.spec.outputFileName = os.path.join( os.getcwd() , SpecBase.genOutputFileName(self._methodName) )
            else:
                if self.spec.outputFileName[-1] == "\\" or self.spec.outputFileName[-1] == "/":
                    self.spec.outputFileName = os.path.join ( os.path.abspath( self.spec.outputFileName ) , SpecBase.genOutputFileName(self._methodName) )

            # ParaMonte variables

            if   self.spec.sampleSize                               is not None: nameList += SpecBase.sampleSize                            (self.spec.sampleSize                         )
            if   self.spec.randomSeed                               is not None: nameList += SpecBase.randomSeed                            (self.spec.randomSeed                         )
            if   self.spec.description                              is not None: nameList += SpecBase.description                           (self.spec.description                        )
            if   self.spec.outputFileName                           is not None: nameList += SpecBase.outputFileName                        (self.spec.outputFileName                     )
            if   self.spec.outputDelimiter                          is not None: nameList += SpecBase.outputDelimiter                       (self.spec.outputDelimiter                    )
            if   self.spec.chainFileFormat                          is not None: nameList += SpecBase.chainFileFormat                       (self.spec.chainFileFormat                    )
            if   self.spec.variableNameList                         is not None: nameList += SpecBase.variableNameList                      (self.spec.variableNameList                   )
            if   self.spec.restartFileFormat                        is not None: nameList += SpecBase.restartFileFormat                     (self.spec.restartFileFormat                  )
            if   self.spec.outputColumnWidth                        is not None: nameList += SpecBase.outputColumnWidth                     (self.spec.outputColumnWidth                  )
            if   self.spec.overwriteRequested                       is not None: nameList += SpecBase.overwriteRequested                    (self.spec.overwriteRequested                 )
            if   self.spec.outputRealPrecision                      is not None: nameList += SpecBase.outputRealPrecision                   (self.spec.outputRealPrecision                )
            if   self.spec.silentModeRequested                      is not None: nameList += SpecBase.silentModeRequested                   (self.spec.silentModeRequested                )
            if   self.spec.domainLowerLimitVec                      is not None: nameList += SpecBase.domainLowerLimitVec                   (self.spec.domainLowerLimitVec                )
            if   self.spec.domainUpperLimitVec                      is not None: nameList += SpecBase.domainUpperLimitVec                   (self.spec.domainUpperLimitVec                )
            if   self.spec.parallelizationModel                     is not None: nameList += SpecBase.parallelizationModel                  (self.spec.parallelizationModel               )
            if   self.spec.progressReportPeriod                     is not None: nameList += SpecBase.progressReportPeriod                  (self.spec.progressReportPeriod               )
            if   self.spec.targetAcceptanceRate                     is not None: nameList += SpecBase.targetAcceptanceRate                  (self.spec.targetAcceptanceRate               )
            if   self.spec.mpiFinalizeRequested                     is not None: nameList += SpecBase.mpiFinalizeRequested                  (self.spec.mpiFinalizeRequested               )
            if   self.spec.maxNumDomainCheckToWarn                  is not None: nameList += SpecBase.maxNumDomainCheckToWarn               (self.spec.maxNumDomainCheckToWarn            )
            if   self.spec.maxNumDomainCheckToStop                  is not None: nameList += SpecBase.maxNumDomainCheckToStop               (self.spec.maxNumDomainCheckToStop            )

            if self._method.isParaDRAM:

                # ParaMCMC variables

                if   self.spec.chainSize                            is not None: nameList += SpecMCMC.chainSize                             (self.spec.chainSize                          )
                if   self.spec.scaleFactor                          is not None: nameList += SpecMCMC.scaleFactor                           (self.spec.scaleFactor                        )
                if   self.spec.startPointVec                        is not None: nameList += SpecMCMC.startPointVec                         (self.spec.startPointVec                      )
                if   self.spec.proposalModel                        is not None: nameList += SpecMCMC.proposalModel                         (self.spec.proposalModel                      )
                if   self.spec.proposalStartCovMat                  is not None: nameList += SpecMCMC.proposalStartCovMat                   (self.spec.proposalStartCovMat                )
                if   self.spec.proposalStartCorMat                  is not None: nameList += SpecMCMC.proposalStartCorMat                   (self.spec.proposalStartCorMat                )
                if   self.spec.proposalStartStdVec                  is not None: nameList += SpecMCMC.proposalStartStdVec                   (self.spec.proposalStartStdVec                )
                if   self.spec.sampleRefinementCount                is not None: nameList += SpecMCMC.sampleRefinementCount                 (self.spec.sampleRefinementCount              )
                if   self.spec.sampleRefinementMethod               is not None: nameList += SpecMCMC.sampleRefinementMethod                (self.spec.sampleRefinementMethod             )
                if   self.spec.randomStartPointRequested            is not None: nameList += SpecMCMC.randomStartPointRequested             (self.spec.randomStartPointRequested          )
                if   self.spec.randomStartPointDomainLowerLimitVec  is not None: nameList += SpecMCMC.randomStartPointDomainLowerLimitVec   (self.spec.randomStartPointDomainLowerLimitVec)
                if   self.spec.randomStartPointDomainUpperLimitVec  is not None: nameList += SpecMCMC.randomStartPointDomainUpperLimitVec   (self.spec.randomStartPointDomainUpperLimitVec)

                # ParaDRAM variables

                if   self.spec.adaptiveUpdateCount                  is not None: nameList += SpecDRAM.adaptiveUpdateCount                   (self.spec.adaptiveUpdateCount                )
                if   self.spec.adaptiveUpdatePeriod                 is not None: nameList += SpecDRAM.adaptiveUpdatePeriod                  (self.spec.adaptiveUpdatePeriod               )
                if   self.spec.greedyAdaptationCount                is not None: nameList += SpecDRAM.greedyAdaptationCount                 (self.spec.greedyAdaptationCount              )
                if   self.spec.delayedRejectionCount                is not None: nameList += SpecDRAM.delayedRejectionCount                 (self.spec.delayedRejectionCount              )
                if   self.spec.burninAdaptationMeasure              is not None: nameList += SpecDRAM.burninAdaptationMeasure               (self.spec.burninAdaptationMeasure            )
                if   self.spec.delayedRejectionScaleFactorVec       is not None: nameList += SpecDRAM.delayedRejectionScaleFactorVec        (self.spec.delayedRejectionScaleFactorVec     )

            nameList = "&" + self._methodName + " " + nameList + SpecBase.interfaceType() + SpecBase.systemInfoFilePath(pm.platform.systemInfoFilePath) + "/"

            ############################################################################################################################
            #### end namelist generation from arguments
            ############################################################################################################################

            inputFileVec_pntr = nameList.encode("utf-8")                # create byte-object from the internal input file

        else:

            if not self.mpiEnabled:
                pm.warn ( msg   = "Input namelist file is given by the user. " + newline
                                + "All simulation specifications will be read from the input file."
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 1
                        )

            inputFileVec_pntr = inputFile.encode("utf-8")                       # create byte-object from the external input file

        inputFileLen = len(inputFileVec_pntr) # byte-object length
        #inputFileLen_pntr = ct.byref( ct.c_size_t( len(inputFileVec_pntr) ) ) # pointer to byte-object length
        inputFileVec_pntr = ct.c_char_p( inputFileVec_pntr )                   # pointer to byte-object

        return inputFileVec_pntr, inputFileLen #_pntr

    ################################################################################################################################
    #### _setFileToRead()
    ################################################################################################################################

    def _setFileToRead(self, file, fileType, fileSuffix):
        if self.spec.outputFileName is None:
            file = os.getcwd()
            if self.reportEnabled:
                pm.warn ( msg   = "The ``file`` is neither given as input to ``read" + fileType.capitalize() + "()``" + newline
                                + "nor set as a simulation specification of the " + self._methodName + " object. " + newline
                                + "This information is essential, otherwise how could the output files be found?" + newline
                                + "All that is needed is the unique name (including path) of the simulation name " + newline
                                + "shared among its output files or simply, the path to the specific " + fileSuffix + newline
                                + " file to be read. For now, the " + self._methodName + " sampler will search " + newline
                                + "the current working directory for simulation output files that match the " + newline
                                + "filename pattern of " + fileSuffix + " files."
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 1
                        )
        else:
            file = self.spec.outputFileName
        return file

    ################################################################################################################################
    #### _setDelimiterToRead()
    ################################################################################################################################

    def _setDelimiterToRead(self, delimiter, fileType, fileSuffix):
        if self.spec.outputDelimiter is None:
            delimiter   = ","
            if self.reportEnabled:
                pm.warn ( msg   = "The ``delimiter`` is neither given as input to ``read" + fileType.capitalize() + "()``" + newline
                                + "nor set as a simulation specification of the " + self._methodName + " object. " + newline
                                + "This information is essential, otherwise how could the output files be parsed?" + newline
                                + "For now, the " + self._methodName + " sampler will assume a comma-separated " + newline
                                + "file format for the contents of the " + fileSuffix + " file(s) to be parsed."
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 1
                        )
        else:
            delimiter = self.spec.outputDelimiter
        return delimiter

    ################################################################################################################################
    #### readTabular()
    ################################################################################################################################

    def _readTabular( self
                    , file          : str
                    , fileType      : str
                    , delimiter     : str
                    , parseContents : bool
                    , renabled      : bool
                    ) -> tp.List[TabularFileContents] :
        """

        Read the contents of the file(s) whose path is given by the input argument ``file``.
        This function is not to be directly accessible to and callable by the users
        of the ParaMonte library.

            **Parameters**

                file

                    A string representing the path to the tabular file with
                    the default value of None.

                    The path only needs to uniquely identify the simulation
                    to which the tabular file belongs. For example, specifying
                    ``"./mydir/mysim"`` as input will lead to a search for a 
                    file that begins with ``"mysim"`` and ends with the tabular file 
                    name's prefix, such as, ``"_sample.txt"``, inside the directory 
                    ``"./mydir/"``. If there are multiple files with such name, 
                    then all of them will be read and returned as a list.

                    If this input argument is not provided by the user, the
                    value of the object attribute ``outputFileName`` will be 
                    used instead. At least one of the two mentioned routes 
                    must provide the path to the tabular file otherwise,
                    this method will break by calling ``sys.exit()``.

                fileType

                    A string containing the type of the file to be parsed.
                    Current options include but are not limited to:
                    ``sample``, ``chain``, ``markovChain``, ``progress``

                delimiter

                    An input string representing the delimiter used in the
                    output tabular file. If it is not provided as input argument,
                    the value of the corresponding object attribute outputDelimiter
                    will be used instead. If none of the two are available,
                    the default comma delimiter "," will be assumed and used.

                parseContents

                    If set to True, the contents of the file will be parsed and
                    stored in a component of the object named ``contents``.
                    The default value is ``True``.

                renabled

                    If set to ``False``, the contents of the file(s) will be 
                    stored as a list in a (new) component of the object with a
                    name that ends with the prefix ``List``. Otherwise, ``None`` 
                    will be the return value of the method. If set to ``True``, 
                    the reverse will done. The default value is ``False``.

            **Returns**

                List

                    A Python list of ``TabularFileContents`` objects, each 
                    of which corresponds to the contents of a unique restart 
                    file. The contents of each object is dependent on the 
                    type of the file that has been parsed.

        """

        if fileType=="sample":
            fileSuffix = "sample"
        elif fileType=="chain" or fileType=="markovChain":
            fileSuffix = "chain"
        elif fileType=="progress":
            fileSuffix = "progress"
        else:
            fileSuffix = None
            pm.abort( msg   = "Internal error occurred. The input fileType is not recognized." + newline
                            + "Please report this error at:" + newline
                            + newline
                            + "    " + pm.website.github.issues.url
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        if file is None: file = self._setFileToRead(file, fileType, fileSuffix)
        if delimiter is None: delimiter = self._setDelimiterToRead(delimiter, fileType, fileSuffix)

        FileList = pm.utils.getFileList(file, fileSuffix, self._methodName, self.reportEnabled)

        tabularContentsList = []
        for file in FileList:

            file = os.path.abspath(file)
            if self.reportEnabled:
                pm.note ( msg = "processing " + fileSuffix + " file: " + file
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 0
                        )

            tabularContents = TabularFileContents   ( file = file
                                                    , fileType = fileType
                                                    , delimiter = delimiter
                                                    , parseContents = parseContents
                                                    , reportEnabled = self.reportEnabled
                                                    , methodName = self._methodName
                                                    )
            tabularContentsList.append(tabularContents)

        outputListName = fileType + "List"
        if renabled:
            outputListFullName = outputListName
            msg=( "The processed " + fileType + " files are now stored in the output variable as a " + newline
                + "Python list. For example, to access the contents of the first (or the only) " + newline
                + fileType + " file stored in an output variable named " + outputListFullName + ", try:"
                )
        else:
            outputListFullName = self._objectName + "." + outputListName
            setattr(self, outputListName, tabularContentsList)
            msg=( "The processed " + fileType + " files are now stored in the newly-created" + newline
                + "component `" + outputListName + "` of the " + self._methodName + " object as a Python list." + newline
                + "For example, to access the contents of the first (or the only) " + fileType + " file, try:"
                )

        if self.reportEnabled:

            if fileSuffix=="progress":
                specials = ""
            else:
                specials =  ( "    " + outputListFullName + "[0].plot.line3()        # to make 3D line plots." + newline
                            + "    " + outputListFullName + "[0].plot.scatter3()     # to make 3D scatter plots." + newline
                            + "    " + outputListFullName + "[0].plot.lineScatter3() # to make 3D line-scatter plots." + newline
                            + "    " + outputListFullName + "[0].plot.contour()      # to make fast 2D kernel density plots." + newline
                            + "    " + outputListFullName + "[0].plot.contourf()     # to make fast 2D kernel density filled contour plots." + newline
                            + "    " + outputListFullName + "[0].plot.contour3()     # to make fast 3D kernel density contour plots." + newline
                            + "    " + outputListFullName + "[0].plot.histplot()     # to make seaborn 1D distribution plots." + newline
                            + "    " + outputListFullName + "[0].plot.kdeplot1()     # to make seaborn 1D kernel density plots." + newline
                            + "    " + outputListFullName + "[0].plot.kdeplot2()     # to make seaborn 2D kernel density plots." + newline
                            + "    " + outputListFullName + "[0].plot.grid()         # to make GridPlot" + newline
                            )

            pm.note ( msg   = msg + newline
                            + newline
                            + "    " + outputListFullName + "[0].df" + newline
                            + newline
                            + "where you will have to replace `pmpd` with your " + self._methodName + " instance name." + newline
                            + "To access the plotting tools, try:" + newline
                            + newline
                            + "    " + outputListFullName + "[0].plot.<PRESS TAB TO SEE THE LIST OF PLOTS>" + newline
                            + newline
                            + "For example," + newline
                            + newline
                            + "    " + outputListFullName + "[0].plot.line()         # to make 2D line plots." + newline
                            + "    " + outputListFullName + "[0].plot.scatter()      # to make 2D scatter plots." + newline
                            + "    " + outputListFullName + "[0].plot.lineScatter()  # to make 2D line-scatter plots." + newline
                            + specials
                            + newline
                            + "To plot or inspect the variable autocorrelations or the correlation/covariance matrices, try:" + newline
                            + newline
                            + "    " + outputListFullName + "[0].stats.<PRESS TAB TO SEE THE LIST OF COMPONENTS>" + newline
                            + newline
                            + "For more information and examples on the usage, visit:" + newline
                            + newline
                            + "    " + pm.website.home.url
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )
        if renabled:
            return tabularContentsList
        else:
            return None

    ################################################################################################################################
    #### readSample
    ################################################################################################################################

    def readSample  ( self
                    , file          : tp.Optional[str]  = None
                    , delimiter     : tp.Optional[str]  = None
                    , parseContents : tp.Optional[bool] = True
                    , renabled      : tp.Optional[bool] = False
                    ) -> tp.List[TabularFileContents] :
        """

        Return a list of the contents of a set of ParaDRAM output
        sample files whose names contain the user-provided input file.
        This method is to be only used for postprocessing of the output
        sample file(s) of an already finished ParaDRAM simulation.
        It is not meant to be called by all processes in parallel mode,
        although it is possible.

            **Parameters**

                file (optional)

                    A string representing the path to the output file with
                    the default value of None.
                    The path only needs to uniquely identify the simulation
                    to which the output file belongs. For example, specifying
                    ``"./mydir/mysim"`` as input will lead to a search for a 
                    file that begins with ``"mysim"`` and ends with ``"_sample.txt"``
                    inside the directory ``"./mydir/"``. If there are multiple
                    files with such name, then all of them will be read
                    and returned as a list.
                    If this input argument is not provided by the user, the
                    value of the object attribute outputFileName
                    will be used instead. At least one of the two mentioned
                    routes must provide the path to the output file otherwise,
                    this method will break by calling ``sys.exit()``.

                delimiter (optional)

                    Optional input string representing the delimiter used in the
                    output output file. If it is not provided as input argument,
                    the value of the corresponding output object's attribute 
                    ``outputDelimiter`` will be used instead. If none of the two 
                    are available, the default comma delimiter ``","`` will be 
                    assumed and used. 

                parseContents (optional)

                    If set to True, the contents of the file will be parsed and
                    stored in a component of the object named ``contents``.
                    The default value is ``True``.

                renabled (optional)

                    If set to ``False``, the contents of the file(s) will be 
                    stored as a list in a (new) component of the ParaDRAM object 
                    named ``sampleList`` and ``None`` will be the return value 
                    of the method. If set to ``True``, the reverse will done.
                    The default value is ``False``.

            **Returns**

                sampleList (optional)

                    A Python list of ``TabularFileContents`` objects, each 
                    of which corresponds to the contents of a unique restart 
                    file. Each object has the following components:

                        file

                            The full absolute path to the output file.

                        delimiter

                            The delimiter used in the output file.

                        ndim

                            The number of dimensions of the domain of the objective 
                            function from which the output has been drawn.

                        count

                            The number of sampled points in the output file.

                        plot

                            A structure containing the graphics tools for the 
                            visualization of the contents of the file.

                        df

                            The contents of the output file in the form of
                            a pandas-library DataFrame (hence called ``df``).

                        contents

                            corresponding to each column in the progress file, a property
                            with the same name as the column header is also created
                            for the object which contains the data stored in that column
                            of the progress file. These properties are all stored in the 
                            attribute ``contents``.

                If ``renabled = True``, the list of objects will be returned as the
                return value of the method. Otherwise, the list will be stored in a
                component of the ParaDRAM object named ``sampleList``.

        """

        return self._readTabular( file = file
                                , fileType = "sample"
                                , delimiter = delimiter
                                , parseContents = parseContents
                                , renabled = renabled
                                )

    ################################################################################################################################
    #### readChain
    ################################################################################################################################

    def readChain   ( self
                    , file          : tp.Optional[str] = None
                    , delimiter     : tp.Optional[str] = None
                    , parseContents : tp.Optional[bool] = True
                    , renabled      : tp.Optional[bool] = False
                    ) -> tp.List[TabularFileContents] :
        """

        Return a list of the contents of a set of ParaDRAM output
        chain files whose names begin the user-provided input file.
        This method is to be only used for postprocessing of the output
        chain file(s) of an already finished ParaDRAM simulation.
        It is not meant to be called by all processes in parallel mode,
        although it is possible.

            **Parameters**

                file (optional)

                    A string representing the path to the output file with
                    the default value of None.
                    The path only needs to uniquely identify the simulation
                    to which the output file belongs. For example, specifying
                    ``"./mydir/mysim"`` as input will lead to a search for a file
                    that begins with ``"mysim"`` and ends with ``"_chain.txt"``
                    inside the directory ``"./mydir/"``. If there are multiple
                    files with such name, then all of them will be read
                    and returned as a list.
                    If this input argument is not provided by the user, the
                    value of the object attribute outputFileName
                    will be used instead. At least one of the two mentioned
                    routes must provide the path to the output file otherwise,
                    this method will break by calling ``sys.exit()``.

                delimiter (optional)

                    Optional input string representing the delimiter used in the
                    output output file. If it is not provided as input argument,
                    the value of the corresponding output object's attribute 
                    ``outputDelimiter`` will be used instead. If none of the two 
                    are available, the default comma delimiter ``","`` will be 
                    assumed and used. 

                parseContents (optional)

                    If set to True, the contents of the file will be parsed and
                    stored in a component of the object named ``contents``.
                    The default value is ``True``.

                renabled (optional)

                    If set to ``False``, the contents of the file(s) will be 
                    stored as a list in a (new) component of the ParaDRAM object 
                    named ``chainList`` and ``None`` will be the return value 
                    of the method. If set to ``True``, the reverse will done.
                    The default value is ``False``.

            **Returns**

                chainList (optional)

                    A Python list of ``TabularFileContents`` objects, each 
                    of which corresponds to the contents of a unique restart 
                    file. Each object has the following components:

                        file

                            The full absolute path to the output file.

                        delimiter

                            The delimiter used in the output file.

                        ndim

                            The number of dimensions of the domain of the objective 
                            function from which the output has been drawn.

                        count

                            The number of sampled points in the output file.

                        plot

                            A structure containing the graphics tools for the 
                            visualization of the contents of the file.

                        df

                            The contents of the output file in the form of
                            a pandas-library DataFrame (hence called ``df``).

                        contents

                            corresponding to each column in the progress file, a property
                            with the same name as the column header is also created
                            for the object which contains the data stored in that column
                            of the progress file. These properties are all stored in the 
                            attribute ``contents``.

                If ``renabled = True``, the list of objects will be returned as the
                return value of the method. Otherwise, the list will be stored in a
                component of the ParaDRAM object named ``sampleList``.

        """

        return self._readTabular( file = file
                                , fileType = "chain"
                                , delimiter = delimiter
                                , parseContents = parseContents
                                , renabled = renabled
                                )

    ################################################################################################################################
    #### readProgress
    ################################################################################################################################

    def readProgress( self
                    , file          : tp.Optional[str] = None
                    , delimiter     : tp.Optional[str] = None
                    , parseContents : tp.Optional[bool] = True
                    , renabled      : tp.Optional[bool] = False
                    ) -> tp.List[TabularFileContents] :
        """

        Return a list of the contents of a set of ParaMonte output
        progress files whose names begin the user-provided input file.
        This method is to be only used for postprocessing of the output
        progress file(s) of an already finished ParaMonte simulation.
        It is not meant to be called by all processes in parallel mode,
        although it is possible.

            **Parameters**

                file (optional)

                    A string representing the path to the output file with
                    the default value of ``None``.

                    The path only needs to uniquely identify the simulation
                    to which the output file belongs. For example, specifying
                    ``"./mydir/mysim"`` as input will lead to a search for a file
                    that begins with ``"mysim"`` and ends with ``"_progress.txt"``
                    inside the directory ``"./mydir/"``. If there are multiple
                    files with such name, then all of them will be read
                    and returned as a list.

                    If this input argument is not provided by the user, the
                    value of the object attribute outputFileName
                    will be used instead. At least one of the two mentioned
                    routes must provide the path to the progress file otherwise,
                    this method will break by calling ``sys.exit()``.

                delimiter (optional)

                    Optional input string representing the delimiter used in the
                    output progress file. If it is not provided as input argument,
                    the value of the corresponding object attribute outputDelimiter
                    will be used instead. If none of the two are available,
                    the default comma delimiter ``","`` will be assumed and used.

                parseContents (optional)

                    If set to True, the contents of the file will be parsed and
                    stored in a component of the object named ``contents``.
                    The default value is ``True``.

                renabled (optional)

                    If set to False, the contents of the file(s) will be stored
                    as a list in a (new) component of the sampler object named
                    ``progressList`` and ``None`` will be the return value of the
                    method. If set to ``True``, the reverse will be done.
                    The default value is ``False``.

            **Returns**

                A list of objects, each of which has the following properties:

                    file

                        The full absolute path to the file.

                    delimiter

                        The delimiter used in the file.

                    ncol

                        The number of columns of the file.

                    plot

                        A structure containing the graphics tools for the 
                        visualization of the contents of the file.

                    df

                        the contents of the progress file in the form of
                        a pandas-library DataFrame (hence called ``df``).

                    contents

                        corresponding to each column in the progress file, a property
                        with the same name as the column header is also created
                        for the object which contains the data stored in that column
                        of the progress file. These properties are all stored in the 
                        attribute ``contents``.

                If ``renabled = True``, the list of objects will be returned as the
                return value of the method. Otherwise, the list will be stored in a
                component of the sampler object named ``progressList``.

        """

        return self._readTabular( file = file
                                , fileType = "progress"
                                , delimiter = delimiter
                                , parseContents = parseContents
                                , renabled = renabled
                                )

    ################################################################################################################################
    #### readRestart
    ################################################################################################################################

    def readRestart ( self
                    , file          : tp.Optional[str] = None
                    , renabled      : tp.Optional[bool] = False
                    ) -> tp.List[RestartFileContents] :
        """

        Return a list of the contents of a set of the simulation(s) 
        output restart files whose names begin the user-provided input 
        file prefix, or as specified by the input simulation specification 
        ``SAMPLER.spec.outputFileName``, where SAMPLER can be an instance 
        of any one of the ParaMonte's sampler classes, such as ``ParaDRAM()``.

            **NOTE**

            Only restart output files in **ASCII format** can be read 
            via this method. The binary restart files are NOT meant to 
            be parsed via this method. To request for ASCII restart 
            output files in simulations, set the input simulation 
            specification

            .. code-block:: python

                SAMPLER.spec.restartFileFormat = "ascii",

            where ``SAMPLER`` can be an instance of any one of the 
            ParaMonte's sampler classes, such as ``ParaDRAM()``.

            **WARNING**
            
            Avoid using this routine for very large long simulations.
            Reading the full restart file of a large-scale simulation 
            problem can be extremely memory-intensive.

            **WARNING**

            This method is to be only used for post-processing of the 
            output restart file(s) of an already finished simulation. 
            It is NOT meant to be called by all processes in parallel 
            mode, although it is possible.

            **Parameters**

                file (optional)

                    A string representing the path to the restart file with the
                    default value of []. The path only needs to uniquely identify
                    the name of the simulation to which the restart file belongs.
                    For example, specifying ``"./mydir/mysim"`` as input will lead 
                    to a search for a file that begins with ``"mysim"`` and ends with
                    ``"_restart.txt"`` inside the directory ``"./mydir/"``.
                    If there are multiple files with such name, then all of them
                    will be read and returned as a list.
                    If this input argument is not provided by the user, the
                    value of the object's ``spec`` attribute ``outputFileName``
                    will be used instead.

                    **WARNING**

                    At least one of the two mentioned routes must
                    provide the path to the restart file. Otherwise,
                    this method will abort the program.

                    Example usage:

                        .. code-block:: python

                            pmpd.readRestart("./out/test_run_")

                        or,

                        .. code-block:: python

                            pmpd.spec.outputFileName = "./out/test_run_"
                            pmpd.readRestart()

                   Both of the above examples are equivalent.
                   The latter is recommended as it is less confusing.

                renabled (optional)

                    If set to ``False``, the contents of the file(s) will be 
                    stored as a list in a (new) component of the sampler object 
                    named ``restartList`` and ``None`` will be the return value 
                    of the method. If set to ``True``, the reverse will done.
                    The default value is ``False``.

            **Returns**

                restartList (optional)

                    A Python list of ``RestartFileContents`` objects, each 
                    of which corresponds to the contents of a unique restart 
                    file. Each object has the following components:

                        file

                            The full absolute path to the restart file.

                        ndim

                            The number of dimensions of the domain of the objective 
                            function for which the restart file was generated.

                        count

                            The number of restart writes to the file.

                        plot

                            A structure containing the graphics tools for the 
                            visualization of the contents of the file.

                        df

                            The contents of the restart file in the form of a
                            ``panda``'s dataframe (``df`` stands for **DataFrame**).

                        contents

                            A structure whose components contain the information 
                            retrieved about each of the entities in the file.

                        propNameList

                            A list of entities names parsed from the restart file.

                If no output argument is provided, a ``restartList`` property will be
                added to the parent sampler-object to which the method ``readRestart()``
                belongs.

        """

        fileType = "restart"
        fileSuffix = "restart"

        if file is None: file = self._setFileToRead(file, fileType, fileSuffix)

        FileList = pm.utils.getFileList(file, fileSuffix, self._methodName, self.reportEnabled)

        restartContentsList = []
        for file in FileList:

            file = os.path.abspath(file)
            if self.reportEnabled:
                pm.note ( msg = "processing " + fileSuffix + " file: " + file
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 0
                        )

            restartContents = RestartFileContents   ( file = file
                                                    , reportEnabled = self.reportEnabled
                                                    , methodName = self._methodName
                                                    )
            restartContentsList.append(restartContents)

        outputListName = fileType + "List"
        if renabled:
            outputListFullName = outputListName
            msg = "The processed " + fileType + " files are now stored in the output variable as a\n" \
                + "Python list. For example, to access the contents of the first (or the only) \n" \
                + fileType + " file stored in an output variable named " + outputListFullName + ", try:"
        else:
            outputListFullName = self._objectName + "." + outputListName
            setattr(self, outputListName, restartContentsList)
            msg = "The processed " + fileType + " files are now stored in the newly-created\n" \
                + "component `" + outputListName + "` of the " + self._methodName + " object as a Python list.\n" \
                + "For example, to access the contents of the first (or the only) " + fileType + " file, try:"

        if self.reportEnabled:
            pm.note ( msg   = msg + newline
                            + newline
                            + "    " + outputListFullName + "[0].contents" + newline
                            + newline
                            + "where you will have to replace `pmpd` with your " + self._methodName + " instance name." + newline
                            + "To access the plotting tools, try:" + newline
                            + newline
                            + "    " + outputListFullName + "[0].plot.<PRESS TAB TO SEE THE LIST OF PLOTS>" + newline
                            + newline
                            + "For example," + newline
                            + newline
                            + "    " + outputListFullName + "[0].plot.line()        # to make bivariate line plots." + newline
                            + "    " + outputListFullName + "[0].plot.scatter()     # to make bivariate scatter plots." + newline
                            + "    " + outputListFullName + "[0].plot.lineScatter() # to make bivariate line-scatter plots." + newline
                            + "    " + outputListFullName + "[0].plot.covmat2()     # to make proposal covariance evolution 2D plots." + newline
                            + "    " + outputListFullName + "[0].plot.covmat3()     # to make proposal covariance evolution 3D plots." + newline
                            + "    " + outputListFullName + "[0].plot.cormat2()     # to make proposal correlation evolution 2D plots." + newline
                            + "    " + outputListFullName + "[0].plot.cormat3()     # to make proposal correlation evolution 3D plots." + newline
                            + newline
                            + "For more information and examples on the usage, visit:" + newline
                            + newline
                            + "    " + pm.website.home.url
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )
        if renabled:
            return restartContentsList
        else:
            return None

    ################################################################################################################################
    #### readReport
    ################################################################################################################################

    def readReport  ( self
                    , file          : tp.Optional[str] = None
                    , renabled      : tp.Optional[bool] = False
                    ) -> tp.List[ReportFileContents] :
        """

        Return a list of the contents of a set of the simulation(s) output
        report files whose names begin the user-provided input file prefix, 
        or as specified by the input simulation specification 
        ``SAMPLER.spec.outputFileName``, where SAMPLER can be 
        an instance of any one of the ParaMonte's sampler 
        classes, such as ``ParaDRAM()``.

        **NOTE**

        This method is to be only used for post-processing of the output
        report file(s) of an already finished simulation. It is NOT meant
        to be called by all processes in parallel mode, although it is
        possible.

            **Parameters**

                file (optional)

                    A string representing the path to the report file with the
                    default value of []. The path only needs to uniquely identify
                    the name of the simulation to which the report file belongs.
                    For example, specifying ``"./mydir/mysim"`` as input will lead to
                    a search for a file that begins with ``"mysim"`` and ends with
                    ``"_report.txt"`` inside the directory ``"./mydir/"``.
                    If there are multiple files with such name, then all of them
                    will be read and returned as a list.
                    If this input argument is not provided by the user, the
                    value of the object's ``spec`` attribute ``outputFileName``
                    will be used instead.

                    **WARNING**

                    At least one of the two mentioned routes must
                    provide the path to the report file. Otherwise,
                    this method will abort the program.

                    Example usage:

                        .. code-block:: python

                            pmpd.readReport("./out/test_run_")

                        or,

                        .. code-block:: python

                            pmpd.spec.outputFileName = "./out/test_run_"
                            pmpd.readReport()

                   Both of the above examples are equivalent.
                   The latter is recommended as it is less confusing.

                renabled (optional)

                    If set to ``False``, the contents of the file(s) will be 
                    stored as a list in a (new) component of the object with a
                    name that ends with the prefix ``List``. Otherwise, ``None`` 
                    will be the return value of the method. If set to ``True``, 
                    the reverse will done. The default value is ``False``.

            **Returns**

                reportList (optional)

                    A Python list of ``ReportFileContents`` objects, each of
                    which corresponds to the contents of a unique report file.
                    Each object may have a dynamic list of the different sections
                    of the output report file
                    Each object may have the following components:

                        file

                            The full absolute path to the report file.

                        contents

                            The contents of the file in its entirely as a string.

                If no output argument is provided, a ``reportList`` property will be
                added to the parent sampler object to which the method ``readReport()``
                belongs.

        """

        fileType = "report"
        fileSuffix = "report"

        if file is None: file = self._setFileToRead(file, fileType, fileSuffix)

        FileList = pm.utils.getFileList(file, fileSuffix, self._methodName, self.reportEnabled)

        reportContentsList = []
        for file in FileList:

            file = os.path.abspath(file)
            if self.reportEnabled:
                pm.note ( msg = "processing " + fileSuffix + " file: " + file
                        , methodName = self._methodName
                        , marginTop = 1
                        , marginBot = 0
                        )

            reportContents = ReportFileContents ( file = file
                                                , reportEnabled = self.reportEnabled
                                                , methodName = self._methodName
                                                )
            reportContentsList.append(reportContents)

        outputListName = fileType + "List"
        if renabled:
            outputListFullName = outputListName
            msg = "The processed " + fileType + " files are now stored in the output variable as a\n" \
                + "Python list. For example, to access the contents of the first (or the only) \n" \
                + fileType + " file stored in an output variable named " + outputListFullName + ", try:"
        else:
            outputListFullName = self._objectName + "." + outputListName
            setattr(self, outputListName, reportContentsList)
            msg = "The processed " + fileType + " files are now stored in the newly-created\n" \
                + "component `" + outputListName + "` of the " + self._methodName + " object as a Python list.\n" \
                + "For example, to access the entire contents of the first (or the only) " + fileType + " file, try:"

        if self.reportEnabled:
            pm.note ( msg   = msg + newline
                            + newline
                            + "    " + outputListFullName + "[0].contents.print()" + newline
                            + newline
                            + "where you will have to replace `pmpd` with your " + self._methodName + " instance name." + newline
                            + "To access the simulation statistics and information, examine the contents of the" + newline
                            + "components of the following structures:" + newline
                            + newline
                            + "    " + outputListFullName + "[0].contents.print()  # to print the contents of the report file." + newline
                            + "    " + outputListFullName + "[0].setup             # to get information about the simulation setup." + newline
                            + "    " + outputListFullName + "[0].stats.time        # to get the timing information of the simulation." + newline
                            + "    " + outputListFullName + "[0].stats.chain       # to get the statistics of the simulation output sample." + newline
                            + "    " + outputListFullName + "[0].stats.numFuncCall # to get information about the number of function calls." + newline
                            + "    " + outputListFullName + "[0].stats.parallelism # to get information about the simulation parallelism." + newline
                            + "    " + outputListFullName + "[0].spec              # to get the simulation specification in the report file." + newline
                            + newline
                            + "For more information and examples on the usage, visit:" + newline
                            + newline
                            + "    " + pm.website.home.url
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )
        if renabled:
            return reportContentsList
        else:
            return None

    ################################################################################################################################
    #### helpme
    ################################################################################################################################

    def helpme  ( self
                , topic : tp.Optional[ str ] = None
                ):
        """

        Prints help on the input object.

            **Parameters**

                topic

                    A string value that is the name of a component of the current 
                    sample object for which help is needed. For example:  

                    Example usage:

                        .. code-block:: python

                            pm.helpme("helpme")

            **Returns**

                None

        """

        usage   = ("    Usage:" + newline
                + newline
                + "        import paramonte as pm " + newline
                + "        pm.helpme()      # to get help on paramonte module. " + newline
                + "        pm.helpme(topic) # to get help on topic. " + newline
                + newline
                + "    where `topic` in the above can be the name of any " + newline
                + "    component of the current sampler object."
                )

        try:
            doc = eval("self."+topic+".__doc__")
            print(doc + "\n")
        except:
            print(self.__doc__)
            print("\nHere is the information on the parent class:\n")
            print(ParaMonteSampler.__doc__)

        if topic.lower()=="helpme": pm.note( msg = usage, methodName = "helpme()", marginTop = 0, marginBot = 1)

        return None

    ################################################################################################################################
