%   +pm - This is the ParaMonte MATLAB parallel Monte Carlo and machine learning library.
%
%   What is ParaMonte?
%   ==================
%
%   ParaMonte is a serial/parallel library of Monte Carlo (**MC**) optimization, sampling, and integration routines
%   for exploring arbitrary-dimensional mathematical density functions, in particular, the posterior distributions of
%   Bayesian models in data science, machine learning, and scientific inference, as well as Machine Leaning (**ML**)
%   algorithms for scientific inference.
%
%   ParaMonte has been designed with the goal of unifying the
%
%       **automation** (of MC and ML simulations),
%       **user-friendliness** (of the library),
%       **accessibility** (from multiple programming environments),
%       **high-performance** (at runtime), and
%       **scalability** (across many parallel processors).
%
%   For more information on the installation, usage, and examples, visit:
%
%       https://www.cdslab.org/paramonte
%
%   To get started with the library, simply create an instance of the @paramonte
%   class to access the routines and functionalities available in the library.
%   Before doing so, ensure to add the library path to your MATLAB session,
%   as done below enclosed between the two comment lines,
%
%       pmlibRootDir = './'; % Change this to the directory containing the ParaMonte MATLAB package.
%       addpath(pmlibRootDir);
%       pm.version
%       doc pm
%
%   Modules
%   -------
%
%       The MATLAB equivalent of modules is called **package**.
%       Similar to modules in other languages, MATLAB packages can be imported to the MATLAB environment
%       or simply be used as namespace. For example,
%
%           pm.lib.version()
%
%       The ParaMonte MATLAB library currently contains the following modules:
%
%           pm.array - module for array manipulation tasks.
%           pm.introspection - module for MATLAB type introspection.
%           pm.io - module for IO tasks.
%           pm.lib - module for ParaMonte library information retrieval.
%           pm.math - module for mathematical tasks.
%           pm.matlab - module for MATLAB-specific tasks.
%           pm.os - module for Operating System-specific tasks.
%           pm.sampling - module for stochastic sampling of mathematical objective functions.
%           pm.sort - module for sorting tasks.
%           pm.str - module for string manipulation tasks.
%           pm.sys - module for system-level operations.
%           pm.timing - module for timing tasks.
%           pm.vis - module for visualization tasks.
%           pm.web - module for world-wide-web tasks.
%
%       Additionally, the package contains the following folder that are meant to be used internally by the library.
%
%           pm/auxil - Folder contains auxiliary information and functionalities required by the library.
%           pm/lib - Folder contains the ParaMonte MATLAB lower-level shared libraries.
%
%   Naming conventions
%   ------------------
%
%       The camelCase naming style is used throughout the entire ParaMonte MATLAB library, across
%       all programming languages including C/C++/Fortran/Julia/MATLAB/Python/R.
%
%       All simulation specifications start with a lowercase letter, including
%       scalar/vector/matrix int, float, string, or Boolean variables.
%
%       Occasionally, variable names representing a vector of values may be suffixed with "Vec".
%
%       Occasionally, variable names representing a matrix of values may be suffixed with "Mat".
%
%       Occasionally, variable names representing a list of values may be suffixed with "List".
%
%       Frequently, functions or class methods begin with a lowercase verb.
%
%       Significant attempt has been made to name all Boolean variables like a logical
%       proposition that is evaluated to either true or false, set by the user.
%
%   Tips
%   ----
%
%       When using the ParaMonte MATLAB library functionalities, particularly ParaMonte samplers in parallel,
%       it would be best to close any such aggressive software/applications as Dropbox, ZoneAlarm, ...
%       that interfere with your ParaMonte MATLAB library output files, potentially causing the tasks to
%       fail and crash before successful completion. These situations happen only scarcely.
%
%       On Windows systems, when restarting an old interrupted ParaDRAM simulation,
%       ensure your MATLAB session is also restarted before the simulation restart.
%       This may be needed as Windows sometimes locks access to some or all of the
%       simulation output files.
%
%       To unset an already-set input simulation specification, simply set the
%       simulation attribute to empty double `[]` or re-instantiate the object.
%
%   LICENSE
%   -------
%
%       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%