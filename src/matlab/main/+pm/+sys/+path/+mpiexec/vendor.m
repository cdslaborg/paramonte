function name = vendor(path)
    %
    %   Return the a MATLAB string containing the MPI library vendor
    %   name corresponding to the input ``mpiexec`` system path.
    %
    %   Parameters
    %   ----------
    %
    %       path
    %
    %           The input scalar MATLAB string,
    %           containing the path to the ``mpiexec``
    %           binary whose vendor is to be determined.
    %           (**optional**,  default = ``pm.sys.path.mpiexec.which()``)
    %
    %   Returns
    %   -------
    %
    %       name
    %
    %           The output MATLAB string containing the MPI library vendor name ALL in lower-case.
    %           Possible values are:
    %
    %               -   ``"Intel"``     :   representing the Intel MPI library.
    %               -   ``"MPICH"``     :   representing the MPICH MPI library.
    %               -   ``"OpenMPI"``   :   representing the OpenMPI library.
    %
    %           If the vendor name cannot be identified, the output will be empty ``[]``.
    %
    %   Interface
    %   ---------
    %
    %       name = pm.sys.path.mpiexec.vendor()
    %       name = pm.sys.path.mpiexec.vendor(path)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    name = "";
    if 0 == nargin
        path = pm.sys.path.mpiexec.which();
        if  path == ""
            return;
        end
    else
        path = string(path);
    end
    if contains(path, ".openmpi")
        name = "OpenMPI";
    elseif contains(path, ".mpich")
        name = "MPICH";
    else
        % Retrieve mpiexec version.
        if isunix()
            path = strrep(path, """", "\""");
        end
        [failed, version] = system("""" + path + """" + " --version");
        versionLower = lower(version);
        if ~failed
            if contains(version, "Intel")
                name = "Intel";
            elseif contains(versionLower, "openmpi") || contains(versionLower, "open-mpi") || contains(versionLower, "openrte")
                name = "OpenMPI";
            elseif contains(versionLower, "mpich") || contains(versionLower, "hydra")
                % MPICH is the most elusive, so we check it as the last possibility.
                name = "MPICH";
            end
        else
            %   On Windows (and perhaps other OS), calling ``system("mpiexec --version")``
            %   can fail because ``mpiexec`` can refer to the MATLAB copy of mpiexec binary.
            %       Unknown option: -version
            %       Error: no executable specified
            %       Unable to parse the mpiexec command arguments.
            %   To bypass this failure, we check if ``cpuinfo`` exists in the same path as ``mpiexec``.
            %   If so, then the mpiexec binary belongs to Intel MPI library.
            [dirname, filename, extname] = fileparts(path);
            if filename == "mpiexec"
                cipath = fullfile(dirname, "cpuinfo" + extname);
                if isfile(cipath)
                    name = "Intel";
                elseif ispc() && extname ~= ".exe"
                    cipath = fullfile(dirname, "cpuinfo.exe");
                    if isfile(cipath)
                        name = "Intel";
                    end
                end
            end
        end
    end
end