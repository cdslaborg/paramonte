%>  \brief
%>  Return the MPI image count and current image ID (e.g., MPI rank + 1)
%>  and the MPI library name as it appears in the ParaMonte library path.<br>
%>
%>  \details
%>  This function attempts to detect the invocation of the current
%>  MATLAB session via an ``mpiexec`` binary MPI launcher instance.<br>
%>
%>  Output values interpretation
%>  ----------------------------
%>
%>  <ol>
%>      <li>    If the output ``mpiname`` is set to an empty string ``""``,
%>              meaning that ``strcmp(mpiname, "")`` yields ``true``,
%>              then the algorithm has failed, possibly implying
%>              that the current system does not support any MPI
%>              parallelism supported by the Paraonte library.
%>      <li>    If the output ``mpiname`` is set to an non-empty string ``""``,
%>              which must be the name of the MPI library vendor used within ParaMonte, then,<br>
%>              then an instance of the corresponding MPI library exists on the current system.<br>
%>              If so,<br>
%>              <ol>
%>                  <li>    A value of ``0`` for either of the output values ``nproc`` and ``rankp1``
%>                          implies the algorithm failed to infer the corresponding parameter.<br>
%>                  <li>    A non-zero positive value for either or both of the output values ``nproc``
%>                          and ``rankp1`` implies that the ``mpiexec`` MPI launcher binary of the
%>                          corresponding MPI library **might have been used** to launch MATLAB
%>                          and its  might have been used to launch MATLAB.<br>
%>                          <ol>
%>                              <li>    A non-zero positive value other than ``1`` for the output ``nproc``
%>                                      guarantees that the ``mpiexec`` launcher has been used to launch MATLAB
%>                                      with the specified ``nproc`` number of MPI images (processes).<br>
%>                          </ol>
%>              </ol>
%>  </ol>
%>
%>  \param[in]  vendor  :   The input scalar MATLAB string, containing the MPI
%>                          library vendor supported by the ParaMonte library.<br>
%>                          Possible values are:<br>
%>                          <ol>
%>                              <li>    ``OpenMPI``, representing the OpenMPI library.
%>                              <li>    ``MPICH``, representing the MPICH MPI library.
%>                              <li>    ``Intel``, representing the Intel MPI library.
%>                              <li>    ``any``, representing any available MPI library.<br>
%>                          </ol>
%>                          or any other informal name returned by the ParaMonte
%>                          MATLAB function [pm.lib.mpi.name()](@ref name).<br>
%>                          Note that **all values are case-insensitive**.<br>
%>                          (**optional, default = ``"any"``.)
%>
%>  \return
%>  ``mpiname``         :   The output scalar MATLAB string containing the MPI
%>                          library name as it appears in the ParaMonte library path.<br>
%>                          An output value of ``""`` implies the algorithm failed to detect
%>                          a ParaMonte-compatible MPI library invoked or installed on the system.<br>
%>  ``bindir``          :   The output scalar MATLAB string containing the directory of the binary within
%>                          the ParaMonte library that succeeded in detecting the local MPI installation.<br>
%>                          This is the directory that contains the ParaMonte MATLAB shared libraries.<br>
%>                          An output value of ``""`` implies the algorithm failed to detect
%>                          a ParaMonte-compatible MPI library invoked or installed on the system.<br>
%>  ``nproc``           :   The output scalar MATLAB non-negative whole number containing
%>                          the number of MPI images in the current invocation of the
%>                          ``mpiexec`` binary MPI launcher.<br>
%>                          An output value of ``0`` implies the algorithm failed to detect
%>                          a ParaMonte-compatible MPI library or failed to infer the image count.<br>
%>  ``rankp1``          :   The output scalar MATLAB non-negative whole number containing
%>                          the image ID (e.g., MPI rank + 1) of the current MPI image in
%>                          the current invocation of the ``mpiexec`` binary MPI launcher.<br>
%>                          An output value of ``0`` implies the algorithm failed to detect
%>                          a ParaMonte-compatible MPI library or failed to infer the image count.<br>
%>                          Note that the image ID always starts at ``1``, unlike the MPI rank.<br>
%>                          The argument ``rankp1`` stands for ``rank + 1``.<br>
%>
%>  \interface{detect}
%>  \code{.m}
%>
%>      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect();
%>      [mpiname, bindir, nproc, rankp1] = pm.lib.mpi.detect(vendor);
%>
%>  \endcode
%>
%>  \warning
%>  This routine must be called before the MPI library is finalized by any other application.<br>
%>  Example applications that finalize the MPI library upon return includes the ParaMonte MPI-parallel samplers.<br>
%>  However, each application offers an option to prevent the default finalization of the MPI library.<br>
%>  For this reason, the results of this function are saved within the function for subsequent usage after the first call.<br>
%>
%>  \note
%>  This function internally uses a ``persistent`` MATLAB ``struct``
%>  to store the results for the default input value for ``vendor``
%>  to expedite future inquiries using this function.<br>
%>  Specifying an input argument value will always
%>
%>  \example{detect}
%>  \include{lineno} example/lib/mpi/detect/main.m
%>  \output{detect}
%>  \include{lineno} example/lib/mpi/detect/main.out.m
%>
%>  \todo
%>  \phigh
%>  An optional ``config`` input argument can be added to this function to allow
%>  to narrow down the list of ParaMonte library paths to search for MPI's presence.<br>
%>
%>  \final{detect}
%>
%>  \author
%>  \AmirShahmoradi, 2:27 AM Monday, November 11, 2024, Dallas, TX<br>
function [mpiname, bindir, nproc, rankp1] = detect(vendor)

    persistent output;

    if  nargin < 1
        vendor = [];
    end

    if  pm.array.len(vendor) == 0
        vendor = "any";
    end

    if  isempty(output) || nargin > 0

        %%%%
        %%%% Get a list of MPI nproc/rank binaries in the current ParaMonte library instance.
        %%%%

        if ~strcmpi(vendor, "any") && ~strcmpi(vendor, "all")
            mpidirname = "**" + pm.lib.mpi.name(vendor);
        else
            mpidirname = "";
        end

        bindirs = pm.sys.path.glob(pm.lib.path.lib() + mpidirname + "**xgetImageStatMPI*");

        %%%%
        %%%% Test the instances for MPI presence.
        %%%%

        output = struct();
        output.nproc = 0;
        output.rankp1 = 0;
        output.bindir = "";
        output.mpiname = "";
        mpinames = ["mpi", "impi", "mpich", "openmpi"];

        for iname = 1 : numel(bindirs)

            %%%%
            %%%% Execute the MPI image status binary.
            %%%%

            [status, cmdout] = system(bindirs{iname});
            cmdout = string(cmdout);

            if  status == 0

                %%%%
                %%%% Find the MPI name in the current ParaMonte library instance.
                %%%%

                [basedir, ~, ~] = fileparts(bindirs{iname});
                [basedir, ~, ~] = fileparts(basedir);
                [~, dirname, ~] = fileparts(string(basedir));
                dirname = string(dirname);

                for choice = mpinames
                    if  strcmp(dirname, choice)
                        output.mpiname = dirname;
                        output.bindir = bindirs(iname);
                        break;
                    end
                end

                %if  strcmp(output.mpiname, "")
                %    warning ( newline ...
                %            + "Unrecognized MPI library name detected among the ParaMonte library paths:" + newline ...
                %            + newline ...
                %            + pm.io.tab + string(output.mpiname) + newline ...
                %            + newline ...
                %            + "The corresponding ParaMonte library path is: " + newline ...
                %            + newline ...
                %            + pm.io.tab + bindirs(iname) + newline ...
                %            + newline ...
                %            + "The expected MPI library names are: " + newline ...
                %            + newline ...
                %            + pm.io.tab + strjoin(mpinames, ", ") + newline ...
                %            + newline ...
                %            );
                %else
                %    return;
                %end

                if  contains(cmdout, "imageCountRank")

                    %%%%
                    %%%% Search for the keyword `imageCountRank`.
                    %%%%

                    tokens = strsplit(cmdout, "imageCountRank");
                    output.rankp1 = str2double(tokens(3));
                    output.nproc = str2double(tokens(2));
                    break;

                end

            end

        end

    end

    mpiname = output.mpiname;
    bindir = output.bindir;
    rankp1 = output.rankp1;
    nproc = output.nproc;

    if  nargin > 0
        output = [];
    end

end