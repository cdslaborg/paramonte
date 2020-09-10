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

import numpy as np
import typing as tp

from _ParaMonteSampler import ParaMonteSampler
from _TabularFileContents import TabularFileContents
import _paramonte as pm

newline = pm.newline

####################################################################################################################################
#### ParaDRAM class
####################################################################################################################################

class ParaDRAM(ParaMonteSampler):
    """

    This is the **ParaDRAM** class to generate instances of **serial** and **parallel**
    **Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo** 
    sampler class of the ParaMonte library. The ``ParaDRAM`` class is a 
    child of the ``ParaMonteSampler`` class.

    All ParaDRAM class attributes (input arguments to the ParaDRAM constructor)
    are optional and all attributes can be also set after a ParaDRAM instance
    is returned by the constructor.

    Once you set the optional attributes to your desired values,
    call the ParaDRAM sampler via the object's method ``runSampler()``.

    .. _example-serial-usage:

    **Example serial usage**

    Copy and paste the following code enclosed between the
    two comment lines in your python/ipython/jupyter session
    (ensure the indentations of the pasted lines comply with Python rules):

    .. code-block:: python
        :linenos:

        ##################################
        import paramonte as pm
        import numpy as np
        def getLogFunc(point):
            # return the log of the standard multivariate
            # Normal density function with ndim dimensions
            return -0.5 * np.sum( np.double( point )**2 )
        pmpd = pm.ParaDRAM()
        pmpd.runSampler ( ndim = 4 # assume 4-dimensional objective function
                        , getLogFunc = getLogFunc   # the objective function
                        )
        ##################################

    where,

        ndim

            represents the number of dimensions of the domain of
            the user's objective function ``getLogFunc(point)`` and,

        getLogFunc(point)

            represents the user's objective function to be sampled,
            which must take a single input argument ``point`` of type
            numpy-float64 array of length ``ndim`` and must return the
            natural logarithm of the objective function.

    .. _example-parallel-usage:

    **Example parallel usage**

    Copy and paste the following code enclosed between the
    two comment lines in your python/ipython/jupyter session
    (ensure the indentations of the pasted lines comply with Python rules):

    .. code-block:: python
        :linenos:

        ##################################
        with open("main.py", "w") as file:
            file.write  ('''
        import paramonte as pm
        import numpy as np
        def getLogFunc(point):
            # return the log of the standard multivariate
            # Normal density function with ndim dimensions
            return -0.5 * np.sum( np.double( point )**2 )
        pmpd = pm.ParaDRAM()
        pmpd.mpiEnabled = True
        pmpd.runSampler ( ndim = 4 # assume 4-dimensional objective function
                        , getLogFunc = getLogFunc   # the objective function
                        )
        ''')
        ##################################

    where,

        ndim

            represents the number of dimensions of the domain of
            the user's objective function ``getLogFunc(point)`` and,

        getLogFunc(point)

            represents the user's objective function that is to be sampled.
            This function must take a single input argument ``point`` of type
            numpy-float64 array of length ndim and must return the natural
            logarithm of the objective function.

        mpiEnabled

            is a logical (boolean) indicator that, if ``True``, will
            cause the ParaDRAM simulation to run in parallel
            on the requested number of processors.
            The default value is ``False``.

    The above will generate a Parallel-ParaDRAM-simulation Python script in the
    current working directory of Python. Note the only difference between the
    serial and parallel simulation scripts: the extra line ``pmpd.mpiEnabled = True``
    which tell the ParaMonte library to run the simulation in parallel.

    Assuming that you already have an MPI runtime library installed on your
    system (see below), you can now execute this Python script file ``main.py``
    in parallel in two ways:

    1.  from inside ipython or jupyter, type the following,

        .. code-block:: bash

            !mpiexec -n 3 python main.py

    2.  outside of Python environment,
        from within a Bash shell (on Linux or Mac) or,
        from within an Anaconda command prompt on Windows,
        type the following,

        .. code-block:: bash

            mpiexec -n 3 python main.py

    **Note:**

    On Windows platform, if you are using the Intel MPI library,
    we recommend that you also specify the extra flag -localonly,

    .. code-block:: bash

        mpiexec -localonly -n 3 python main.py

    This will cause the simulations to run in parallel only on a single node,
    but more importantly, it will also prevent the use of Hydra service and
    the requirement for its registration. If you are not on a Windows cluster,
    (e.g., you are using your personal device), then we highly recommend
    specifying this flag.


    In all cases in the above, the script ``main.py`` will run on 3 processors.
    Feel free to change the number of processors to any number desired. But do
    not request more than the available number of physical cores on your system.

    **Tips on parallel usage**

    For up-to-date detailed instructions on how to run simulations in parallel visit:

        https://www.cdslab.org/paramonte

    You can also use the following commands on the Python command-line,

    .. code-block:: python
        :linenos:

        ##################################
        import paramonte as pm
        pm.verify() # verify the existence of parallel simulation prerequisites
        ##################################

    to obtain specific information on how to run a parallel simulation,
    in particular, in relation to your current installation of ParaMonte.
    In general, for parallel simulations:

    0.  Ensure you need and will get a speedup by running the ParaDRAM sampler
        in parallel. Typically, if a single evaluation of the objective function
        takes much longer than a few milliseconds, your simulation may then
        benefit from the parallel run.

    1.  Ensure you have an MPI library installed, preferably, the Intel MPI
        runtime libraries. An MPI library should be automatically installed
        on your system with ParaMonte. If needed, you can download the Intel
        MPI library from their website and install it.

    2.  Ensure the ParaDRAM object property ``mpiEnabled`` is ``True``
        (the default is ``False``).

    3.  Before running the parallel simulation, in particular, on Windows systems,
        you may need to define the necessary MPI environmental variables on your system.
        To get information on how to define the variables, use the paramonte module's
        function, ``verify()``, as described in the above.

    4.  Call your main Python code from a Python-aware mpiexec-aware command-line via,

        .. code-block:: bash

            mpi_launcher -n num_process python name_of_yor_python_code.py

        where,

        1.  "mpi_launcher" is the name of the MPI launcher
            of the MPI runtime library that you have installed.
            For example, the Intel MPI library's launcher is named mpiexec,
            also recognized by Microsoft, MPICH, and OpenMPI.
            Note that on supercomputers, the MPI launcher is usually
            something other than ``mpiexec``, for example:
            ``ibrun``, ``mpirun``, ...

        2.  "num_process" represents the number of cores
            on which you want to run the program. Replace this
            with the an integer number, like, 3 (meaning 3 cores).

            Do not assign more processes than the available number of
            physical cores on your device/cluster. Assigning more cores
            than physically available on your system will only slow down
            your simulation.

    Once the above script is saved in the file ``main.py``, open a Python-aware and
    MPI-aware command prompt to run the simulation in parallel via the MPI launcher,

    .. code-block:: bash

        mpiexec -n 3 python main.py

    This will execute the Python script ``main.py`` on three processes (images).
    Keep in mind that on Windows systems you may need to define MPI environmental
    variables before a parallel simulation, as described in the above.

    **ParaDRAM Class Attributes**

    See also:

        https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/

    All input specifications (attributes) of a ParaDRAM simulation are optional.
    However, it is recommended that you provide as much information as possible
    about the specific ParaDRAM simulation and the objective function to be sampled
    via ParaDRAM simulation specifications.

    The ParaDRAM simulation specifications have lengthy comprehensive descriptions
    that appear in full in the output report file of every ParaDRAM simulation.

    The best way to learn about individual ParaDRAM simulation attributes
    is to a run a minimal serial simulation with the following Python script,

    .. code-block:: python
        :linenos:

        ##################################
        from paramonte import ParaDRAM
        pmpd = ParaDRAM()
        pmpd.spec.outputFileName = "./test"
        def getLogFunc(point): return -sum(point**2)
        pmpd.runSampler( ndim = 1, getLogFunc = getLogFunc )
        ##################################

    Running this code will generate a set of simulation output files (in the current
    working directory of Python) that begin with the prefix ``test_process_1``. Among
    these, the file ``test_process_1_report.txt`` contains the full description of all
    input specifications of the ParaDRAM simulation as well as other information
    about the simulation results and statistics.

    **Parameters**

        None. The simulation specifications can be set once an object is instantiated.
        All simulation specification descriptions are collectively available at:

            https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/

        Note that this is the new interface. The previous ParaDRAM class interface
        used to optionally take all simulation specifications as input. However,
        overtime, this approach has become more of liability than any potential
        benefit. All simulation specifications have to be now to be set solely
        after a ParaDRAM object is instantiated, instead of setting the
        specifications via the ParaDRAM class constructor.

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
            If it is set to ``True``, it will cause the ParaDRAM simulation
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

            Specifying an input file will cause the ParaDRAM sampler
            to ignore all other simulation specifications set by the 
            user via sampler instance's `spec`-component attributes.

        spec

            A frozen class containing all simulation specifications.
            All simulation attributes are by default set to appropriate
            values at runtime. To override the default simulation
            specifications, set the `spec` attributes to some
            desired values of your choice. For possible values, see:

                https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/

            If you need help on any of the simulation specifications, try
            the supplied ``helpme()`` function in this component, like,

            .. code-block:: python
                :linenos:

                ##################################
                import paramonte as pm
                pmpd = pm.ParaDRAM()          # instantiate a ParaDRAM sampler class
                pmpd.spec.helpme()            # get help on all simulation specification
                pmpd.spec.helpme("chainSize") # get help on "chainSize" specifically
                ##################################

    **Methods**

        See below for information on the methods.

    **Returns**

        Object of class ParaDRAM sampler.

    ---------------------------------------------------------------------------
    """

    def __init__(self):
        """

        The constructor for ParaDRAM class.
        All input parameters are optional and all class
        attributes can be changed after the object construction.

            **Parameters**

                None

        """

        super().__init__(methodName = "ParaDRAM")

        ## ParaMonte specifications
        #
        #self.spec = pm.utils.FrozenClass()
        #
        ## ParaMonte variables
        #self.spec.sampleSize                           = sampleSize
        #self.spec.randomSeed                           = randomSeed
        #self.spec.description                          = description
        #self.spec.outputFileName                       = outputFileName
        #self.spec.outputDelimiter                      = outputDelimiter
        #self.spec.chainFileFormat                      = chainFileFormat
        #self.spec.variableNameList                     = variableNameList
        #self.spec.restartFileFormat                    = restartFileFormat
        #self.spec.outputColumnWidth                    = outputColumnWidth
        #self.spec.outputRealPrecision                  = outputRealPrecision
        #self.spec.silentModeRequested                  = silentModeRequested
        #self.spec.domainLowerLimitVec                  = domainLowerLimitVec
        #self.spec.domainUpperLimitVec                  = domainUpperLimitVec
        #self.spec.parallelizationModel                 = parallelizationModel
        #self.spec.progressReportPeriod                 = progressReportPeriod
        #self.spec.targetAcceptanceRate                 = targetAcceptanceRate
        #self.spec.mpiFinalizeRequested                 = mpiFinalizeRequested
        #self.spec.maxNumDomainCheckToWarn              = maxNumDomainCheckToWarn
        #self.spec.maxNumDomainCheckToStop              = maxNumDomainCheckToStop
        ## ParaMCMC variables
        #self.spec.chainSize                            = chainSize
        #self.spec.scaleFactor                          = scaleFactor
        #self.spec.startPointVec                        = startPointVec
        #self.spec.proposalModel                        = proposalModel
        #self.spec.proposalStartCovMat                  = proposalStartCovMat
        #self.spec.proposalStartCorMat                  = proposalStartCorMat
        #self.spec.proposalStartStdVec                  = proposalStartStdVec
        #self.spec.sampleRefinementCount                = sampleRefinementCount
        #self.spec.sampleRefinementMethod               = sampleRefinementMethod
        #self.spec.randomStartPointRequested            = randomStartPointRequested
        #self.spec.randomStartPointDomainLowerLimitVec  = randomStartPointDomainLowerLimitVec
        #self.spec.randomStartPointDomainUpperLimitVec  = randomStartPointDomainUpperLimitVec
        ## ParaDRAM variables
        #self.spec.adaptiveUpdateCount                  = adaptiveUpdateCount
        #self.spec.adaptiveUpdatePeriod                 = adaptiveUpdatePeriod
        #self.spec.greedyAdaptationCount                = greedyAdaptationCount
        #self.spec.delayedRejectionCount                = delayedRejectionCount
        #self.spec.burninAdaptationMeasure              = burninAdaptationMeasure
        #self.spec.delayedRejectionScaleFactorVec       = delayedRejectionScaleFactorVec
        #
        #self.spec.helpme = SpecDRAM.helpme
        #self.spec._freeze()

    ################################################################################################################################
    #### runSampler
    ################################################################################################################################

    def runSampler  ( self
                    , ndim          : int
                    , getLogFunc    : tp.Callable[[tp.List[float]], float]
                    , inputFile     : tp.Optional[str] = None
                    ) -> None:
        """

        Run ParaDRAM sampler and return nothing.

            **Parameters**

                ndim

                    An integer representing the number of dimensions of the
                    domain of the user's objective function ``getLogFunc(point)``.
                    It must be a positive integer.

                getLogFunc(point)

                    represents the user's objective function to be sampled,
                    which must take a single input argument ``point`` of type
                    numpy-float64 array of length ``ndim`` and must return the
                    natural logarithm of the objective function.

                inputFile (optional)

                    A string input representing the path to an external
                    input namelist of simulation specifications.

                        **WARNING**

                        Use this optional argument with caution and only
                        if you know what you are doing. Specifying this option
                        will cause the sampler to ignore all other simulation
                        specifications set by the user via the ``spec``
                        component of the sampler instance.

            **Returns**

                None

        """

        if not isinstance(ndim,int) or ndim<1:
            pm.abort( msg   = "The input argument ndim must be a positive integer," + newline
                            + "representing the number of dimensions of the domain of" + newline
                            + "the user's objective function getLogFunc()." + newline
                            + "You have entered ndim = " + str(ndim)
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        if not callable(getLogFunc):
            pm.abort( msg   = "The input argument getLogFunc must be a callable function." + newline
                            + "It represents the user's objective function to be sampled," + newline
                            + "which must take a single input argument of type numpy" + newline
                            + "float64 array of length ndim and must return the" + newline
                            + "natural logarithm of the objective function."
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        if inputFile is not None and not isinstance(inputFile,str):
            pm.abort( msg   = "The input argument ``inputFile`` must be of type str." + newline
                            + "It is an optional string input representing the path to" + newline
                            + "an external input namelist of simulation specifications." + newline
                            + "USE THIS OPTIONAL ARGUMENT WITH CAUTION AND" + newline
                            + "ONLY IF YOU KNOW WHAT YOU ARE DOING." + newline
                            + "Specifying this option will cause the sampler to ignore" + newline
                            + "all other simulation specifications set by the user via" + newline
                            + "the ``spec`` component of the sampler instance." + newline
                            + "You have entered inputFile = " + str(inputFile)
                    , methodName = self._methodName
                    , marginTop = 1
                    , marginBot = 1
                    )

        def getLogFunc2arg(ndim,point):
            PointVec = np.array(point[0:ndim])
            return getLogFunc(PointVec)

        self._runSampler( ndim
                        , getLogFunc2arg
                        , inputFile
                        )

    ################################################################################################################################
    #### readMarkovChain
    ################################################################################################################################

    def readMarkovChain ( self
                        , file          : tp.Optional[str] = None
                        , delimiter     : tp.Optional[str] = None
                        , parseContents : tp.Optional[bool] = True
                        , renabled      : tp.Optional[bool] = False
                        ) -> tp.List[TabularFileContents] :
        """

        Return a list of the unweighted verbose (Markov-chain) contents
        of a set of ParaDRAM output chain files, whose names begin the
        user-provided input variable ``file``. This method is to be only
        used for the postprocessing of the output chain file(s) of an
        already finished ParaDRAM simulation. It is not meant to be called
        by all processes in parallel mode, although it is possible.

            **Parameters**

                file (optional)

                    A string representing the path to the chain file with
                    the default value of ``None``.
                    The path only needs to uniquely identify the simulation
                    to which the chain file belongs. For example, specifying
                    ``"./mydir/mysim"`` as input will lead to a search for a file
                    that begins with ``"mysim"`` and ends with ``"_chain.txt"``
                    inside the directory ``"./mydir/"``. If there are multiple
                    files with such name, then all of them will be read
                    and returned as a list.
                    If this input argument is not provided by the user, the
                    value of the object attribute ``outputFileName`` will be
                    used instead. At least one of the two mentioned routes
                    must provide the path to the chain file otherwise,
                    this method will break by calling ``sys.exit()``.

                delimiter (optional)

                    An input string representing the delimiter used in the
                    output chain file. If it is not provided as input argument, the
                    value of the corresponding object attribute ``outputDelimiter``
                    will be used instead. If none of the two are available,
                    the default comma delimiter ``","`` will be assumed and used.

                parseContents (optional)

                    If set to ``True``, the contents of the file will be parsed
                    and stored in a component of the object named ``contents``.
                    The default value is ``True``.

                renabled (optional)

                    If set to False, the contents of the file(s) will be stored
                    as a list in a (new) component of the ParaDRAM object named
                    ``markovChainList`` and ``None`` will be the return value
                    of the method. If set to True, the reverse will done.
                    The default value is ``False``.

            **Returns**

                A list of objects, each of which has the following properties:

                    file

                        The full absolute path to the chain file.

                    delimiter

                        The delimiter used in the chain file.

                    ndim

                        The number of dimensions of the domain of the objective
                        function from which the chain has been drawn.

                    count

                        The number of unique (weighted) points in the chain file.
                        This is essentially the number of rows in the chain file
                        minus one (representing the header line).

                    plot

                        A structure containing the graphics tools for the 
                        visualization of the contents of the file.

                    df

                        The unweighted (Markovian) contents of the chain file in the
                        form of a pandas-library DataFrame (hence called ``df``).

                    contents

                        corresponding to each column in the progress file, a property
                        with the same name as the column header is also created
                        for the object which contains the data stored in that column
                        of the progress file. These properties are all stored in the 
                        attribute ``contents``.

                If ``renabled = True``, the list of objects will be returned as the
                return value of the method. Otherwise, the list will be stored in a
                component of the ParaDRAM object named ``markovChainList``.

        """

        return self._readTabular( file = file
                                , fileType = "markovChain"
                                , delimiter = delimiter
                                , parseContents = parseContents
                                , renabled = renabled
                                )

    ################################################################################################################################

