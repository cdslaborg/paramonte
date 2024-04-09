function names = choices()
    %
    %   Return the ParaMonte-preferred MPI library vendor name(s)
    %   as used in naming the ParaMonte MATLAB shared libraries
    %   in the order of preference on the current platform.
    %
    %   \note
    %
    %       Only the default (first) mpi library name is guaranteed to be
    %       supported in any pre-built distribution of the ParaMonte library.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       names
    %
    %           The output vector of MATLAB strings containing the
    %           the ParaMonte-preferred MPI library vendor names as
    %           used in naming the ParaMonte MATLAB shared libraries.
    %           in the default order of preference.
    %
    %   Interface
    %   ---------
    %
    %       names = pm.lib.mpi.choices()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if ispc()
        names = "impi";
    elseif pm.os.is.lin()
        names = ["impi", "mpich", "openmpi"];
    elseif ismac()
        names = ["openmpi", "mpich"];
    end
end