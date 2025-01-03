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
    bindirs = pm.sys.path.glob(pm.lib.path.lib() + "**pm_sampling*");
    mpinames = ["mpi", "impi", "mpich", "openmpi"];
    names = strings(numel(bindirs), 1);
    iname = 0;
    for ibin = 1 : numel(bindirs)
        [basedir, ~, ~] = fileparts(bindirs(ibin));
        [basedir, ~, ~] = fileparts(basedir);
        [~, dirname, ~] = fileparts(string(basedir)); % This is the library checking mode.
        [~, dirname, ~] = fileparts(string(dirname)); % This is the parallelism name used in the folder naming.
        if  any(strcmp(["mpi", "impi", "mpich", "openmpi"], dirname))
            iname = iname + 1;
            names(iname) = string(dirname);
        elseif ~any(strcmp(["serial", "openmp"], dirname))
            warning ( newline ...
                    + "Unrecognized MPI library name detected among the ParaMonte library paths:" + newline ...
                    + newline ...
                    + pm.io.tab() + string(dirname) + newline ...
                    + newline ...
                    + "The corresponding ParaMonte library path is: " + newline ...
                    + newline ...
                    + pm.io.tab() + bindirs(ibin) + newline ...
                    + newline ...
                    + "The expected MPI library names are: " + newline ...
                    + newline ...
                    + pm.io.tab() + strjoin(mpinames, ", ") + newline ...
                    + newline ...
                    );
        end
    end
    names = unique(names(1 : iname));
    %if  ispc()
    %    names = "impi";
    %elseif pm.os.is.lin()
    %    names = ["impi", "mpich", "openmpi"];
    %elseif ismac()
    %    names = ["openmpi", "mpich"];
    %end
end