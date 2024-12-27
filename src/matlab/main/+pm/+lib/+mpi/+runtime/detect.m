%>  \brief
%>  Return the MPI image count and current image ID (e.g., MPI rank + 1)
%>  and the MPI library name as it appears in the ParaMonte library path
%>  and the corresponding MPI-parallel ParaMonte library path .<br>
%>
%>  \details
%>  This function attempts to detect the invocation of the
%>  current MATLAB session via checking the runtime behavior
%>  of the ``mpiexec`` binary MPI launcher instance.<br>
%>
%>  Output values interpretation
%>  ----------------------------
%>
%>  <ol>
%>      <li>    If the output ``mpiname`` is set to an empty string ``""``,
%>              meaning that ``strcmp(mpiname, "")`` yields ``true``,
%>              then the algorithm has failed, implying that the current
%>              application is likely to have been launched by via an ``mpiexec``
%>              binary whose MPI library is supported by the ParaMonte library.<br>
%>      <li>    If the output ``mpiname`` is set to an non-empty string,
%>              which must be the name of an MPI library vendor used within ParaMonte,<br>
%>              then an instance of the corresponding MPI library likely exists on the current system
%>              and its ``mpiexec`` binary has been used to launch the current application in parallel.<br>
%>              If so,<br>
%>              <ol>
%>                  <li>    A value of ``0`` for either of the output values ``nproc`` and ``rankp1``
%>                          implies the algorithm failed to infer the corresponding parameters.<br>
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
%>      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect();
%>      [mpiname, nproc, rankp1] = pm.lib.mpi.runtime.detect(vendor);
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
%>  \include{lineno} example/lib/mpi/runtime/detect/main.m
%>  \output{detect}
%>  \include{lineno} example/lib/mpi/runtime/detect/main.out.m
%>
%>  \todo
%>  \phigh
%>  The internals of this function can be improved for better efficiency and less redundancy.<br>
%>
%>  \final{detect}
%>
%>  \author
%>  \AmirShahmoradi, 2:27 AM Monday, November 11, 2024, Dallas, TX<br>
function [mpiname, nproc, rankp1] = detect(vendor)

    if  nargin < 1
        vendor = [];
    end

    if  pm.array.len(vendor) == 0
        vendor = "any";
    end

    %%%%
    %%%% Get the name list of available ParaMonte-compatible MPI library names.
    %%%%

    allnames = pm.lib.mpi.choices();
    isall = strcmpi(vendor, "any") || strcmpi(vendor, "all");
    if ~isall
        thisname = pm.lib.mpi.name(vendor);
        isall = strcmpi(thisname, "mpi");
    end
    if  isall
        mpinames = allnames;
    else
        mpinames = thisname;
    end

    %%%%
    %%%% Get the name list of available ParaMonte-compatible MPI library names.
    %%%%

    nproc = 0;
    rankp1 = 0;
    mpiname = "";
    for iname = 1 : numel(allnames)
        thisname = allnames(iname);
        if  any(contains(mpinames, thisname))
            if  strcmp(thisname, "openmpi")
                [nproc, rankp1] = pm.lib.mpi.runtime.ompi();
            elseif strcmp(thisname, "mpich")
                [nproc, rankp1] = pm.lib.mpi.runtime.mmpi();
            elseif strcmp(thisname, "impi")
                [nproc, rankp1] = pm.lib.mpi.runtime.impi();
            end
            if  nproc ~= 0 && rankp1 ~= 0
                mpiname = thisname;
                return;
            end
        end

    end

    % persistent output;
    % output = struct();
    %
    % if  nargin < 1
    %     vendor = [];
    % end
    %
    % if  pm.array.len(vendor) == 0
    %     vendor = "any";
    % end
    %
    % %%%%
    % %%%% Get a list of MPI nproc/rank binaries in the current ParaMonte library instance.
    % %%%%
    %
    % %mpinames = ["mpi", "impi", "mpich", "openmpi"];
    % if ~strcmpi(vendor, "any") && ~strcmpi(vendor, "all")
    %     mpidirpattern = "**" + pm.lib.mpi.name(vendor);
    % else
    %     mpidirpattern = "";
    % end
    %
    % bindirs = pm.sys.path.glob(pm.lib.path.lib() + mpidirpattern + "**xgetImageStatMPI*");
    %
    % %%%%
    % %%%% Loop over the binaries to test them for MPI compatibility.
    % %%%% Beware that some MPI libraries may be compatible but not function correctly when mixed.
    % %%%% Therefore, we will test for the MPI world size larger than ``1`` as an MPI detection.
    % %%%%
    %
    % nprocmax = -1;
    % for iname = 1 : numel(bindirs)
    %
    %     %%%%
    %     %%%% Find the MPI name in the current ParaMonte library instance.
    %     %%%%
    %
    %     [basedir, ~, ~] = fileparts(bindirs{iname});
    %     [basedir, ~, ~] = fileparts(basedir);
    %     [~, dirname, ~] = fileparts(string(basedir));
    %     dirname = string(dirname);
    %
    %     isfresh = true;
    %     if ~isempty(fields(output))
    %         for field = string(fields(output))
    %             isfresh = ~strcmp(dirname, field);
    %             if ~isfresh
    %                 break;
    %             end
    %         end
    %     end
    %
    %     %%%%
    %     %%%% Test the MPI vendor compatibility only if it is a fresh search.
    %     %%%%
    %
    %     if  isfresh
    %
    %         output.(dirname) = struct();
    %         output.(dirname).mpiname = "";
    %         output.(dirname).bindir = "";
    %         output.(dirname).rankp1 = 0;
    %         output.(dirname).nproc = 0;
    %
    %         %%%%
    %         %%%% Execute the MPI image status binary.
    %         %%%%
    %
    %         [status, cmdout] = system(bindirs{iname});
    %         cmdout = string(cmdout);
    %
    %         if  status == 0
    %
    %             %if  strcmp(output.mpiname, "")
    %             %    warning ( newline ...
    %             %            + "Unrecognized MPI library name detected among the ParaMonte library paths:" + newline ...
    %             %            + newline ...
    %             %            + pm.io.tab() + string(output.mpiname) + newline ...
    %             %            + newline ...
    %             %            + "The corresponding ParaMonte library path is: " + newline ...
    %             %            + newline ...
    %             %            + pm.io.tab() + bindirs(iname) + newline ...
    %             %            + newline ...
    %             %            + "The expected MPI library names are: " + newline ...
    %             %            + newline ...
    %             %            + pm.io.tab() + strjoin(mpinames, ", ") + newline ...
    %             %            + newline ...
    %             %            );
    %             %else
    %             %    return;
    %             %end
    %
    %             %%%%
    %             %%%% Search for the keyword `imageCountRank`.
    %             %%%%
    %
    %             if  contains(cmdout, "imageCountRank")
    %
    %                 tokens = strsplit(cmdout, "imageCountRank");
    %                 output.(dirname).nproc = str2double(tokens(2));
    %                 output.(dirname).rankp1 = str2double(tokens(3));
    %                 output.(dirname).bindir = bindirs(iname);
    %                 output.(dirname).mpiname = dirname;
    %                 %break;
    %
    %             end
    %
    %         end
    %
    %     end
    %
    %     %%%%
    %     %%%% Select the MPI library with the highest ``nproc``.
    %     %%%%
    %
    %     if  nprocmax < output.(dirname).nproc
    %         nprocmax = output.(dirname).nproc;
    %         selection = output.(dirname).mpiname;
    %     end
    %
    % end
    %
    % mpiname = output.(selection).mpiname;
    % bindir = output.(selection).bindir;
    % rankp1 = output.(selection).rankp1;
    % nproc = output.(selection).nproc;
    %
    % if  nargin > 0
    %     output = [];
    % end

end