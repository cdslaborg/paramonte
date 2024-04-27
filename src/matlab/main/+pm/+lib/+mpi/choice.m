function name = choice()
    %
    %   Return the ParaMonte-preferred MPI library vendor name
    %   as used in naming the ParaMonte MATLAB shared libraries.
    %   Support for the output MPI library name by this routine is
    %   guaranteed on the current platform.
    %
    %   Parameters
    %   ----------
    %
    %       None
    %
    %   Returns
    %   -------
    %
    %       name
    %
    %           The output scalar MATLAB string containing the
    %           the ParaMonte-preferred MPI library vendor name as
    %           used in naming the ParaMonte MATLAB shared libraries.
    %
    %   Interface
    %   ---------
    %
    %       name = pm.lib.mpi.choice()
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    if ispc() || pm.os.is.lin()
        name = "impi";
    elseif ismac()
        name = "openmpi";
    end
end