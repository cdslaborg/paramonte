%>  \brief
%>  Return the a MATLAB string containing the MPI library vendor
%>  name corresponding to the input ``mpiexec`` system path.<br>
%>
%>  \param[in]  path    :   The input scalar MATLAB string,
%>                          containing the path to the ``mpiexec``
%>                          binary whose vendor is to be determined.<br>
%>                          (**optional**,  default = [pm.sys.path.mpiexec.which()](@ref which))
%>
%>  \return
%>  ``name``            :   The output MATLAB string containing the MPI library vendor name ALL in lower-case.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``"Intel"``     :   representing the Intel MPI library.
%>                              <li>    ``"MPICH"``     :   representing the MPICH MPI library.
%>                              <li>    ``"OpenMPI"``   :   representing the OpenMPI library.
%>                          </ol>
%>                          If the vendor name cannot be identified, the output will be empty ``""``.<br>
%>
%>  \interface{vendor}
%>  \code{.m}
%>
%>      name = pm.sys.path.mpiexec.vendor()
%>      name = pm.sys.path.mpiexec.vendor(path)
%>
%>  \endcode
%>
%>  \see
%>  [pm.lib.mpi.name()](@ref name)<br>
%>  [pm.lib.mpi.vendor()](@ref vendor)<br>
%>  [pm.sys.path.mpiexec.vendor()](@ref vendor)<br>
%>
%>  \example{vendor}
%>  \include{lineno} example/sys/path/mpiexec/vendor/main.m
%>  \output{vendor}
%>  \include{lineno} example/sys/path/mpiexec/vendor/main.out.m
%>
%>  \final{vendor}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:04 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function name = vendor(path)
    name = "";
    if  0 == nargin
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
        failed = pm.os.is.lin() && startsWith(path, "/mnt/"); % A Windows-path application in WSL freezes MATLAB.
        if ~failed
            [failed, version] = system("""" + path + """" + " --version");
            failed = failed ~= 0;
        end
        if ~failed
            versionLower = lower(version);
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