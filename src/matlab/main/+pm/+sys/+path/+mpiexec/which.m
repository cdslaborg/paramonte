%>  \brief
%>  Return the a MATLAB string containing the path to
%>  the first ``mpiexec`` executable binary found in system path.
%>
%>  \details
%>  The output of this function is what is typically returned
%>  by the Unix command ``command -v mpiexec`` or the Windows CMD
%>  command ``where mpiexec``.
%>
%>  \warning
%>  Note that this function intentionally skips any ``mpiexec``
%>  path that is installed by and within the MATLAB binary directory.
%>
%>  \param[in]  vendor  :   The input scalar MATLAB string, containing the MPI
%>                          library vendor that should match the ``mpiexec`` binary.
%>                          Possible values are:<br>
%>                          -   "Intel"     :   representing the Intel MPI library.
%>                          -   "MPICH"     :   representing the MPICH MPI library.
%>                          -   "OpenMPI"   :   representing the OpenMPI library.
%>                          (**optional**,  default = ``""`` representing all possible vendors.)
%>
%>  \return
%>  `path`              :   The output MATLAB string containing the paths to
%>                          the first ``mpiexec`` executable binary found in system path.
%>                          If the ``mpiexec`` does not exist or match the specified ``vendor``,
%>                          the output will be an empty string ``""``.
%>
%>
%>  \interface{which}
%>  \code{.m}
%>
%>      path = pm.sys.path.mpiexec.which()
%>      path = pm.sys.path.mpiexec.which(vendor)
%>
%>  \endcode
%>  \final{which}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 5:06 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function path = which(vendor)
    if 0 < nargin
        vendorLower = lower(vendor);
        if vendorLower == "impi"
            vendorLower = "intel";
        elseif vendorLower == "open-mpi" || vendorLower == "openrte"
            vendorLower = "openmpi";
        end
    else
        vendorLower = "";
    end
    path = "";
    appNameList = ["mpiexec", "mpiexec.openmpi", "mpiexec.mpich"];
    for appName = appNameList
        if ispc()
            [failed, stdout] = system("where " + appName);
        else
            [failed, stdout] = system("command -v " + appName);
        end
        if ~failed
            path = string(strip(stdout, newline));
            %   On Windows, ``where()`` returns a list of all identified paths.
            %   Note that MATLAB has also its own ``mpiexec`` which can show up in ``path``.
            if ispc()
                pathList = split(path, newline);
                for i = 1 : length(pathList)
                    item = pathList(i);
                    if ~contains(item, "MATLAB" + filesep + "R")
                        path = pathList(i);
                    end
                end
            end
            vendorNameLower = lower(pm.sys.path.mpiexec.vendor(path));
            if  vendorNameLower ~= ""
                if ~(contains(appName, vendorLower) || contains(vendorNameLower, vendorLower))
                    path = "";
                end
            else
                path = "";
            end
        end
        if  path ~= ""
            path = string(strip(path, newline));
            break;
        end
    end
end