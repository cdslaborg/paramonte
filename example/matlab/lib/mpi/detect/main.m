cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect()')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect()
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("impi")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("impi")
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("INTEL")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("INTEL")
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("intel")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("intel")
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("ompi")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("ompi")
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("openmpi")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("openmpi")
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("open-mpi")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("open-mpi")
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("mpich")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("mpich")
disp("")

disp('[mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("mpi")')
      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect("mpi")