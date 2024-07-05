%>  \brief
%>  Return a list of scalar MATLAB strings containing the paths to
%>  all detected ``mpiexec`` binaries installed on the system and
%>  available in the environment path variable.<br>
%>
%>  \details
%>  The search strategy is to parse and search the paths in the environmental
%>  ``PATH`` variable of the runtime processor shell and return all ``mpiexec`` paths.<br>
%>  Also, all ``mpiexec`` paths found via [pm.sys.path.mpiexec.which(vendor)](@ref which) are returned.<br>
%>  Additionally, if ``vendor`` is missing or is set to ``"Intel"``, also search the default
%>  installation directories of Intel MPI libraries on all operating systems.<br>
%>  Think of this functionality ``pm.sys.path.mpiexec.find(vendor)`` as a more
%>  comprehensive of what [pm.sys.path.mpiexec.which(vendor)](@ref which) does.<br>
%>
%>  \warning
%>  In Microsoft Windows Subsystem for Linux (WSL) environments, this function
%>  may freeze MATLAB indefinitely, if Windows paths containing MPI installations
%>  exist within the environmental variable ``PATH``.<br>
%>  Such paths typically begin with ``/mnt/`` and are automatically
%>  prepended by WSL to the contents of the ``PATH`` variable.<br>
%>  This problem happens because Windows applications cannot be executed
%>  from within a WSL terminal, freezing the terminal indefinitely.<br>
%>  The best solution is to remove Windows paths from the contents of the
%>  environment variable ``PATH`` of the WSL terminals.<br>
%>  See [this StackOverflow question](https://stackoverflow.com/questions/51336147/how-to-remove-the-win10s-path-from-wsl)
%>  for possible solutions.<br>
%>
%>  \param[in]  vendor  :   The input scalar MATLAB string, containing the MPI
%>                          library vendor that should match the ``mpiexec`` binary.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``Intel``, representing the Intel MPI library.
%>                              <li>    ``MPICH``, representing the MPICH MPI library.
%>                              <li>    ``OpenMPI``, representing the OpenMPI library.
%>                          </ol>
%>                          (**optional**,  default = ``""``)
%>
%>  \return
%>  ``pathList``        :   A list of scalar MATLAB strings containing the paths to
%>                          all detected ``mpiexec`` binaries installed on the system.<br>
%>                          If the ``mpiexec`` is not found or does not match the specified ``vendor``,
%>                          the output will be an empty list ``[]``.<br>
%>
%>  \interface{find}
%>  \code{.m}
%>
%>      pathList = pm.sys.path.mpiexec.find()
%>      pathList = pm.sys.path.mpiexec.find(vendor)
%>
%>  \endcode
%>
%>  \example{find}
%>  \include{lineno} example/sys/path/mpiexec/find/main.m
%>  \output{find}
%>  \include{lineno} example/sys/path/mpiexec/find/main.out.m
%>
%>  \final{find}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:03 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function pathList = find(vendor)

    pathList = [];
    if  0 < nargin
        vendorLower = lower(vendor);
        if  vendorLower == "impi"
            vendorLower = "intel";
        elseif vendorLower == "open-mpi" || vendorLower == "openrte"
            vendorLower = "openmpi";
        end
        vendorList = vendorLower;
    else
        vendorList = ["Intel", "MPICH", "OpenMPI"];
        vendorLower = "";
    end

    %%%%
    %%%% The strategy is to search for any executable in the
    %%%% environmental paths whose name contains `mpiexec`.
    %%%%

    paths = getenv("PATH");
    paths = string(strsplit(paths, pathsep));
    for path = paths
        apps = pm.sys.path.list(path);
        for icell = 1 : length(apps)
            app = apps(icell);
            if  contains(app, "mpiexec") && isfile(app)
                if  0 < nargin
                    name = lower(pm.sys.path.mpiexec.vendor(app));
                    if  name ~= ""
                        if  contains(name, vendorLower)
                            pathList = [pathList, app];
                        end
                    end
                else
                    pathList = [pathList, app];
                end
            end
        end
    end

    %%%%
    %%%% Try more via `which()` and default installation paths.
    %%%%

    for name = vendorList
        path = pm.sys.path.mpiexec.which(name);
        if  path ~= ""
            if  isempty(pathList)
                pathList = [pathList, path];
            elseif sum(strcmp(pathList, path)) == 0
                pathList = [pathList, path];
            end
        end
    end

    if  nargin == 0 || contains(vendorLower, "intel")
        possibilities = [ "C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin\mpiexec.exe" ...
                        , "C:\Program Files\Intel\oneAPI\mpi\latest\bin\mpiexec.exe" ...
                        , "/opt/intel/oneapi/mpi/latest/bin/mpiexec" ...
                        ];
        for path = possibilities
            if  isfile(path)
                if isempty(pathList)
                    pathList = [pathList, path];
                elseif sum(contains(pathList, path)) == 0
                    pathList = [pathList, path];
                end
            end
        end
    end

end