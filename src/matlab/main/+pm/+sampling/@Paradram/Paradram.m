%>  \brief
%>  This is the ParaDRAM class for generating instances of serial and parallel
%>  Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo
%>  sampler of the ParaMonte MATLAB library.
%>
%>  \details
%>  Once you assign the desired simulation specifications to the corresponding
%>  attributes within the component `spec` of an object of class ``pm.sampling.Paradram``,
%>  call the ParaDRAM sampler via the object method [pm.sampling.Paradram.run()](@ref Paradram::run).<br>
%>
%>  While the constructor of this class does not take any input arguments,
%>  all ParaDRAM simulation specifications can be set after creating the object.
%>
%>  \return
%>  `sampler`   :   The output scalar object of class [pm.sampling.Paradram](@ref Paradram).
%>
%>  \interface{Paradram}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>
%>  \endcode
%>
%>  \warning
%>  When using the ParaMonte MATLAB library functionalities, particularly ParaMonte samplers in parallel,
%>  it would be best to close any such aggressive software/applications as **Dropbox**, **ZoneAlarm**, ...
%>  that interfere with your ParaMonte MATLAB library output files, potentially causing the tasks
%>  to fail and crash before successful completion.<br>
%>  These situations scarcely happen.<br>
%>
%>  \note
%>  On Windows systems, when restarting an old interrupted ParaDRAM simulation,
%>  ensure your MATLAB session is also restarted before the simulation restart.<br>
%>  This may be needed as **Windows frequently locks access to some or all simulation output files**.<br>
%>
%>  \note
%>  To unset an already-set input simulation specification, simply set the
%>  simulation attribute to empty double `[]` or re-instantiate an object of class [pm.sampling.Paradram.run()](@ref Paradram).<br>
%>
%>  \see
%>  [ParaDRAM simulation specifications listing](\pmdoc_usage_sampling/paradram/specifications/)<br>
%>  [ParaDRAM simulation restart functionality](\pmdoc_usage_sampling/paradram/restart/)<br>
%>  [ParaDRAM simulation output files](\pmdoc_usage_sampling/paradram/output/)<br>
%>
%>  \example{himmelblau}
%>  \include{lineno} example/sampling/Paradram/himmelblau/main.m
%>  \include{lineno} example/sampling/Paradram/himmelblau/main.out.m
%>  \vis{himmelblau}
%>  \image html example/sampling/Paradram/himmelblau/domain.png width=700
%>  \image html example/sampling/Paradram/himmelblau/traceplot.png width=700
%>  \image html example/sampling/Paradram/himmelblau/proposalAdaptation.png width=700
%>
%>  Example Usage: Serial
%>  ---------------------
%>
%>  First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.
%>
%>  Here is a MATLAB script ``main.m`` for a serial ParaDRAM simulation.<br>
%>  Copy and paste the following code into your MATLAB session:<br>
%>
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>      sampler.run ( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function.
%>                  , 4                 ... ndim:       the number of dimensions of the objective function.
%>                  );
%>      samples = sampler.readSample();
%>      sample = samples{1};
%>      tile = pm.vis.tile.Line(sample.df);
%>      tile.make("coly", sample.sampleLogFuncColIndex + 1 : sample.sampleLogFuncColIndex + 4, "colc", "sampleLogFunc");
%>
%>  \endcode
%>
%>  The mathematical objective function in the above example is a
%>  is a multivariate Normal distribution centered at the origin,
%>  whose natural logarithm is returned by the lambda (``Anonymous``)
%>  function defined as a function handle input to the ParaDRAM sampler.<br>
%>
%>  Running this code will generate a set of simulation output files (in the current working directory of MATLAB).<br>
%>  Among these, the file suffixed with "_report.txt" contains the full description of all input specifications
%>  of the ParaDRAM simulation as well as other information about the simulation results.<br>
%>
%>  Example Usage: Thread-Parallel
%>  ------------------------------
%>
%>  First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.<br>
%>
%>  Threading parallelism is possible as of ParaMonte MATLAB version ``3.0.0``.<br>
%>  However, only ``singleChain`` ParaDRAM simulations are supported.<br>
%>
%>  Here is a MATLAB script ``main.m`` for a thread-parallel ParaDRAM simulation.<br>
%>  Copy and paste the following code and paste into your MATLAB session:<br>
%>
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>      sampler.spec.parallelismNumThread = 0; % use all available threads.
%>      sampler.run ( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function.
%>                  , 4                 ... ndim:       the number of dimensions of the objective function.
%>                  );
%>      samples = sampler.readSample();
%>      sample = samples{1};
%>      pm.vis.tile(sample.contents)
%>
%>  \endcode
%>
%>  The mathematical objective function in the above example is a
%>  is a multivariate Normal distribution centered at the origin,
%>  whose natural logarithm is returned by the lambda (``Anonymous``)
%>  function defined as a function handle input to the ParaDRAM sampler.<br>
%>
%>  Running this code will generate a set of simulation output files (in the current working directory of MATLAB).<br>
%>  Among these, the file suffixed with ``"_report.txt"`` contains the full description of all input specifications
%>  of the ParaDRAM simulation as well as other information about the simulation results.<br>
%>
%>  Specifying ``0`` as the number of threads will lead to using
%>  all available CPU threads for thread-parallel ParaDRAM simulation.<br>
%>
%>  \note
%>  **Benefits of thread-parallelism**<br>
%>  Thread-parallel simulations offer a much more flexible
%>  and easier approach to benefiting from parallelism without
%>  going through the hassle of MPI-parallel simulations.<br>
%>  But they can still potentially offer much faster speed than serial simulations.<br>
%>  The actual speedup depends on a lot of factors.<br>
%>  Moreover, the number of threads is limited to maximum
%>  number of physical cores available on your system.<br>
%>  As such, thread-parallel simulations are not scalable.<br>
%>  If you need scalability, checkout MPI-parallelism below.<br>
%>
%>  Example Usage: MPI-Parallel
%>  ---------------------------
%>
%>  First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.<br>
%>
%>  MPI-parallel simulations can be slightly more cumbersome than thread-parallel simulations
%>  described above because MPI-parallel simulations cannot be performed from within a MATLAB GUI
%>  session and require launching MATLAB via a compatible ``mpiexec`` launcher.<br>
%>
%>  <ol>
%>      <li>    Ensure you need and will get a speedup by running the an MPI-parallel simulation.<br>
%>              Typically, your simulation may then benefit from parallelism only if a single
%>              evaluation of the objective function takes longer than a few milliseconds.<br>
%>
%>      <li>    Ensure the required MPI libraries are installed on your system
%>              (You can skip this step if you know that you already have
%>              a compatible MPI library installed on your system).<br>
%>              On the MATLAB command line type the following,<br>
%>              \code{.m}
%>                  pm.lib.verify();
%>              \endcode
%>              This will verify the existence of a valid MPI library on your system and,
%>              if missing, will guide you to install the MPI library on your system.<br>
%>
%>      <li>    Once the MPI installation is verified, copy and paste the following
%>              code into your MATLAB session:
%>              \code{.m}
%>
%>                  fid = fopen("main_mpi.m", "w");
%>                  sourceCode = ...
%>                  "sampler = pm.sampling.Paradram();" + newline + ...
%>                  "sampler.mpiname = pm.lib.mpi.choice();" + newline + ...
%>                  "sampler.run( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function." + newline + ...
%>                  "           , 4                 ... ndim:       the number of dimensions of the objective function" + newline + ...
%>                  "           );";
%>                  fprintf(fid, "%s\n", sourceCode);
%>                  fclose(fid);
%>
%>              \endcode
%>
%>      <li>    This will generate a ``main_mpi.m`` MATLAB script file in the current working directory of your MATLAB session.<br>
%>              Now, you can execute this MATLAB script file (``main_mpi.m``) in parallel.<br>
%>              To do so, you need to call MATLAB on a command-line, **out of MATLAB GUI**.<br>
%>              <ol>
%>                  <li>    **On Windows**:<br>
%>                          From within command prompt that recognizes both MATLAB and ``mpiexec``,
%>                          ideally, the Intel dedicated command-prompt that is shipped with Intel MPI library,
%>                          type the following,
%>                          \code{.m}
%>
%>                              mpiexec -localonly -n 3 matlab -batch "main_mpi"
%>
%>                          \endcode
%>
%>                          \note
%>                          In the above MPI launcher command for Windows OS,
%>                          we assumed that you would be using the Intel MPI library, hence,
%>                          the reason for the extra flag ``-localonly``.<br>
%>                          This flag runs the parallel code only on one node, but in doing so,
%>                          it avoids the use of Hydra service and its registration.<br>
%>                          If you are not on a Windows cluster, (e.g., you are using your personal device),
%>                          then we recommend specifying this flag.<br>
%>
%>                  <li>    **On macOS/Linux**:<br>
%>                          From within a Bash terminal that recognizes both MATLAB and ``mpiexec``,
%>                          type the following,
%>                          \code{.m}
%>
%>                              mpiexec -n 3 matlab -batch "main_mpi"
%>
%>                          \endcode
%>
%>                          \note
%>                          In both cases in the above, the script ``main_mpi.m`` will run on 3 processors.<br>
%>                          Feel free to change the number of processors to any number desired.<br>
%>                          But do not request more than the available number of physical cores on your system.<br>
%>              </ol>
%>  </ol>
%>
%>  \warning
%>  Do not add postprocessing codes (such as reading and plotting the output samples) in your MPI-parallel code.<br>
%>  There is no point in doing so, since MATLAB will run in ``-batch`` mode for parallel simulations, disabling all plotting capabilities.<br>
%>  Moreover, if you read and postprocess the output files in parallel mode, the task will be done
%>  by all parallel processes, potentially overwriting external IO activities of each other.<br>
%>  Only perform the sampling as described above in MPI-parallel mode.<br>
%>
%>  ParaDRAM Simulation Specifications
%>  ----------------------------------
%>
%>  The ParaDRAM simulation specifications have lengthy comprehensive descriptions
%>  that appear in full in the output report files of every ParaDRAM simulation.<br>
%>
%>  The best way to learn about individual ParaDRAM simulation attributes
%>  is to a run a minimal serial simulation as given in the above.<br>
%>  You can also use the ``sampler.spec.doc()`` method:
%>
%>  \code{.m}
%>      sampler = pm.sampling.Paradram();
%>      sampler.spec.doc();
%>  \endcode
%>
%>  \final{Paradram}
%>
%>  \author
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin
classdef Paradram < pm.sampling.Sampler
    methods(Access = public)
        function self = Paradram()
            self = self@pm.sampling.Sampler("ParaDRAM")
        end
        failed = run(self, getLogFunc, ndim);
        chainMarkovList = readChainMarkov(self, pattern, sep);
    end
end