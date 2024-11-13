cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect()')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect()
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("impi")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("impi")
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("INTEL")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("INTEL")
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("intel")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("intel")
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("ompi")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("ompi")
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("openmpi")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("openmpi")
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("open-mpi")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("open-mpi")
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("mpich")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("mpich")
disp("")

disp('[mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("mpi")')
      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect("mpi")