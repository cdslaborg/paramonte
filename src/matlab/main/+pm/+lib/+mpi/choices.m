%>  \brief
%>  Return the ParaMonte-preferred MPI library vendor name(s)
%>  as used in naming the ParaMonte MATLAB shared libraries
%>  in the order of preference on the current platform.
%>
%>  \return
%>  ``names``   :   The output vector of MATLAB strings containing the
%>                  the ParaMonte-preferred MPI library vendor names as
%>                  used in naming the ParaMonte MATLAB shared libraries.<br>
%>
%>  \interface{choices}
%>  \code{.m}
%>
%>      names = pm.lib.mpi.choices()
%>
%>  \endcode
%>
%>  \warning
%>  Beware that this function only verifies the existence of the
%>  relevant MPI-parallel ParaMonte libraries on the system.<br>
%>  It does not detect the existence of the relevant MPI
%>  runtime library installations on the system.<br>
%>
%>  \example{choices}
%>  \include{lineno} example/lib/mpi/choices/main.m
%>  \output{choices}
%>  \include{lineno} example/lib/mpi/choices/main.out.m
%>
%>  \final{choices}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:18 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function names = choices()
    bindirs = pm.sys.path.glob(pm.lib.path.lib() + "**pm_parallelism*");
    mpinames = ["mpi", "impi", "mpich", "openmpi"];
    names = strings(numel(bindirs), 1);
    for iname = 1 : numel(bindirs)
        [basedir, ~, ~] = fileparts(bindirs(iname));
        [basedir, ~, ~] = fileparts(basedir);
        [~, dirname, ~] = fileparts(string(basedir));
        names(iname) = string(dirname);
        if ~any(strcmp(["mpi", "impi", "mpich", "openmpi"], names(iname)))
            warning ( newline ...
                    + "Unrecognized MPI library name detected among the ParaMonte library paths:" + newline ...
                    + newline ...
                    + pm.io.tab + string(names(iname)) + newline ...
                    + newline ...
                    + "The corresponding ParaMonte library path is: " + newline ...
                    + newline ...
                    + pm.io.tab + bindirs(iname) + newline ...
                    + newline ...
                    + "The expected MPI library names are: " + newline ...
                    + newline ...
                    + pm.io.tab + strjoin(mpinames, ", ") + newline ...
                    + newline ...
                    );
        end
    end
    names = unique(names);
    %if  ispc()
    %    names = "impi";
    %elseif pm.os.is.lin()
    %    names = ["impi", "mpich", "openmpi"];
    %elseif ismac()
    %    names = ["openmpi", "mpich"];
    %end
end