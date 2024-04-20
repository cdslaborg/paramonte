classdef Paradram < pm.sampling.Sampler
    %
    %   This is the ParaDRAM class for generating instances of serial and parallel
    %   Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo
    %   sampler of the ParaMonte MATLAB library.
    %
    %   Once you assign the desired simulation specifications to the corresponding
    %   attributes within the component `spec` of an object of class ``pm.sampling.Paradram``,
    %   call the ParaDRAM sampler via the object method ``pm.sampling.Paradram.run()``.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %           See also the documentation of the class constructor.
    %
    %       \note
    %
    %           All simulation specifications can be set after creating the object.
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
    %       The output scalar object of class ``pm.sampling.Paradram``.
    %
    %   Interface
    %   ---------
    %
    %       sampler = pm.sampling.Paradram();
    %
    %   Example Usage: Serial
    %   ---------------------
    %
    %       First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.
    %
    %       Here is a MATLAB script ``main.m`` for a serial ParaDRAM simulation.
    %       Copy and paste the following code enclosed between the
    %       two comment lines in your MATLAB session:
    %
    %           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           sampler = pm.sampling.Paradram();
    %           sampler.run ( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function.
    %                       , 4                 ... ndim:       the number of dimensions of the objective function.
    %                       );
    %           samples = sampler.readSample();
    %           sample = samples{1};
    %           pm.vis.tile(sample.contents)
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
    %   Example Usage: Thread-Parallel
    %   ------------------------------
    %
    %       First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.
    %
    %       Threading parallelism is possible as of ParaMonte MATLAB version 3.
    %       However, only ``singleChain`` ParaDRAM simulations are supported.
    %
    %       Here is a MATLAB script ``main.m`` for a thread-parallel ParaDRAM simulation.
    %       Copy and paste the following code enclosed between the
    %       two comment lines in your MATLAB session:
    %
    %           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           sampler = pm.sampling.Paradram();
    %           sampler.run ( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function.
    %                       , 4                 ... ndim:       the number of dimensions of the objective function.
    %                       );
    %           samples = sampler.readSample();
    %           sample = samples{1};
    %           pm.vis.tile(sample.contents)
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
    %       Specifying ``0`` as the number of threads will lead to using
    %       all available CPU threads for thread-parallel ParaDRAM simulation.
    %
    %       Benefits of thread-parallelism
    %
    %       Thread-parallel simulations offer a much more flexible
    %       and easier approach to benefiting from parallelism without
    %       going through the hassle of MPI-parallel simulations,
    %       yet, potentially much faster than serial simulations.
    %       The actual speedup depends on a lot of factors.
    %       Moreover, the number of threads is limited to maximum
    %       number of physical cores available on your system.
    %       As such, thread-parallel simulations as not scalable.
    %       If you need scalability, checkout MPI-parallelism below.
    %
    %   Example Usage: MPI-Parallel
    %   ---------------------------
    %
    %       First, ensure the ParaMonte ``+pm`` package (i.e., folder) is available in your MATLAB paths.
    %
    %       MPI-parallel simulations can be slightly more cumbersome than thread-parallel simulations
    %       described above because MPI-parallel simulations cannot be performed from within a MATLAB GUI
    %       session and require launching MATLAB via a compatible ``mpiexec`` launcher.
    %
    %       0.  Ensure you need and will get a speedup by running the an MPI-parallel simulation.
    %           Typically, your simulation may then benefit from parallelism only if a single
    %           evaluation of the objective function takes longer than a few milliseconds.
    %
    %       1.  First, ensure the required MPI libraries are installed on your System:
    %           (You can skip this step if you know that you already have
    %           a compatible MPI library installed on your system).
    %           On the MATLAB command line type the following,
    %
    %               pm.lib.verify();
    %
    %           This will verify the existence of a valid MPI library on your system and,
    %           if missing, will guide you to install the MPI library on your system.
    %
    %       2.  Once the MPI installation is verified,
    %           copy and paste the following code enclosed
    %           between the two comment lines in your MATLAB session:
    %
    %               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               fid = fopen("main_mpi.m", "w");
    %               sourceCode = ...
    %               "sampler = pm.sampling.Paradram();" + newline + ...
    %               "sampler.mpiname = pm.lib.mpi.choice();" + newline + ...
    %               "sampler.run( @(x) - sum(x.^2)  ... getLogFunc: the natural log of the objective function." + newline + ...
    %               "           , 4                 ... ndim:       the number of dimensions of the objective function" + newline + ...
    %               "           );";
    %               fprintf(fid, "%s\n", sourceCode);
    %               fclose(fid);
    %               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       3.  This will generate a ``main_mpi.m`` MATLAB script file in the current
    %           working directory of your MATLAB session. Now, you can execute
    %           this MATLAB script file (``main_mpi.m``) in parallel. To do so,
    %           you need to call MATLAB on a command-line, **out of MATLAB GUI**.
    %
    %               a.  On Windows:
    %
    %                   From within command prompt that recognizes both MATLAB and ``mpiexec``,
    %                   ideally, the Intel dedicated command-prompt that is shipped with Intel MPI library,
    %                   type the following,
    %
    %                       mpiexec -localonly -n 3 matlab -batch "main_mpi"
    %
    %                   \note
    %
    %                       In the above MPI launcher command for Windows OS,
    %                       we assumed that you would be using the Intel MPI library, hence,
    %                       the reason for the extra flag ``-localonly``. This flag runs the parallel
    %                       code only on one node, but in doing so, it avoids the use of Hydra service
    %                       and its registration. If you are not on a Windows cluster, (e.g., you are
    %                       using your personal device), then we recommend specifying this flag.
    %
    %               b.  On macOS/Linux:
    %
    %                   From within a Bash terminal that recognizes both MATLAB and ``mpiexec``,
    %                   type the following,
    %
    %                       mpiexec -n 3 matlab -batch "main_mpi"
    %
    %           \note
    %
    %               In both cases in the above, the script ``main_mpi.m`` will run on 3 processors.
    %               Feel free to change the number of processors to any number desired. But do not
    %               request more than the available number of physical cores on your system.
    %
    %       \warning
    %
    %           Do not add postprocessing codes (such as reading and plotting the output samples)
    %           in your MPI-parallel code. There is no point in doing so, since MATLAB will run in `-batch`
    %           mode for parallel simulations, disabling all plotting capabilities. Moreover, if you read
    %           and postprocess the output files in parallel mode, the task will be done by all parallel
    %           processes, potentially overwriting external IO activities of each other.
    %           Only perform the sampling as described above in MPI-parallel mode.
    %
    %   ParaDRAM Simulation Specifications
    %   ----------------------------------
    %
    %       The ParaDRAM simulation specifications have lengthy comprehensive descriptions
    %       that appear in full in the output report files of every ParaDRAM simulation.
    %
    %       The best way to learn about individual ParaDRAM simulation attributes
    %       is to a run a minimal serial simulation as given in the above.
    %       You can also use the ``sampler.spec.doc()`` method:
    %
    %           sampler = pm.sampling.Paradram();
    %           sampler.spec.doc();
    %
    %       See also: https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    methods(Access = public)
        function self = Paradram()
            %
            %   Construct and return a scalar object of class ``pm.sampling.Paradram``.
            %
            %   Parameters
            %   ----------
            %
            %       None
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Interface
            %   ---------
            %
            %       sampler = pm.sampling.Paradram();
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            self = self@pm.sampling.Sampler("ParaDRAM")
        end
        failed = run(self, getLogFunc, ndim);
        chainMarkovList = readChainMarkov(self, pattern, sep);
    end
end