function str = name(vendor)
    %
    %   Return the MPI library name as used in
    %   naming the ParaMonte MATLAB shared libraries.
    %
    %   Parameters
    %   ----------
    %
    %       vendor
    %
    %           The input scalar MATLAB string, containing the MPI
    %           library vendor supported by the ParaMonte library.
    %           Possible values are:
    %
    %               -   ``Intel``, representing the Intel MPI library.
    %               -   ``MPICH``, representing the MPICH MPI library.
    %               -   ``OpenMPI``, representing the OpenMPI library.
    %
    %           All values are case-insensitive.
    %
    %   Returns
    %   -------
    %
    %       str
    %
    %           The output scalar MATLAB string containing the MPI
    %           library name corresponding to the input MPI ``vendor``.
    %           If the input ``vendor`` is not supported, the default
    %           string ``"mpi"`` is returned as the value ``str``.
    %
    %   Interface
    %   ---------
    %
    %       str = pm.lib.mpi.name(vendor)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    vendorLower = lower(string(vendor));
    for i = length(vendorLower) : -1 : 1
        if strcmp(vendorLower(i), "openmpi") || strcmp(vendorLower(i), "open-mpi") || strcmp(vendorLower(i), "ompi")
            str(i) = "openmpi";
        elseif strcmp(vendorLower(i), "mpich") || strcmp(vendorLower(i), "mmpi")
            str(i) = "mpich";
        elseif strcmp(vendorLower(i), "intel") || strcmp(vendorLower(i), "impi")
            str(i) = "impi";
        else
            str(i) = "mpi";
        end
    end
end