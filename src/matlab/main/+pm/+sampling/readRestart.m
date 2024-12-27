%>  \brief
%>  Return a list of objects of superclass [pm.sampling.FileContentsRestart](@ref FileContentsRestart)
%>  containing the contents of a (set of) ParaMonte simulation output restart file(s) whose paths match
%>  the specified input ``pattern`` or the simulation specification ``sampler.spec.outputFileName``.<br>
%>
%>  \warning
%>  This function is to be only used for post-processing of the output restart file(s) of an already finished simulation.<br>
%>  Although possible, this method is **not** meant to be called by all processes in MPI-parallel simulations.<br>
%>
%>  \param[in]  sampler :   The input object of superclass [pm.sampling.Sampler](@ref Sampler)
%>                          whose type and properties determine the type of the output object(s).<br>
%>                          (**optional**. If empty, it is set to [pm.sampling.Sampler()](@ref Sampler).)
%>  \param[in]  pattern :   The input scalar MATLAB string containing the pattern matching
%>                          the desired restart file(s) whose contents is to be read.<br>
%>                          The specified ``pattern`` only needs to partially identify
%>                          the name of the simulation to which the restart file belongs.<br>
%>                          For example, specifying ``"./mydir/mysim"`` as input will
%>                          lead to a search for file(s) beginning with ``"mysim"`` and
%>                          ending with ``"_restart.txt"`` inside the directory ``"./mydir/"``.<br>
%>                          If there are multiple files matching in the input ``pattern``,
%>                          then all such files will be read and returned as elements of a list.<br>
%>                          If the specified pattern is a valid existing URL, the file will be
%>                          downloaded as a temporary file to the local system, its contents
%>                          shall be parsed and the file will be subsequently removed.<br>
%>                          If the input ``pattern`` is empty, then the method will search
%>                          for any possible candidate files with the appropriate suffix
%>                          in the current working directory.<br>
%>                          (**optional**. If empty, it is set to ``sampler.spec.outputFileName`` or if empty, it is set to ``"./"``.)
%>
%>  \return
%>  ``restartList``     :   The output MATLAB cell array of objects of superclass [pm.sampling.FileContentsRestart](@ref FileContentsRestart),
%>                          each of which corresponds to the contents of a unique ParaMonte sampler restart file.<br>
%>                          <ol>
%>                              <li>    If the input argument ``sampler`` is of type [pm.sampling.Sampler](@ref Sampler),
%>                                      then the output array elements are of type [pm.sampling.FileContentsRestart](@ref FileContentsRestart).<br>
%>                              <li>    If the input argument ``sampler`` is of type [pm.sampling.Paradram](@ref Paradram),
%>                                      then the output array elements are of type [pm.sampling.FileContentsRestartDRAM](@ref FileContentsRestartDRAM).<br>
%>                          </ol>
%>
%>  \interface{readRestart}
%>  \code{.m}
%>
%>      restartList = pm.sampling.readRestart();
%>      restartList = pm.sampling.readRestart([], []);
%>      restartList = pm.sampling.readRestart([], pattern);
%>      restartList = pm.sampling.readRestart(sampler, []);
%>      restartList = pm.sampling.readRestart(sampler, pattern);
%>
%>  \endcode
%>
%>  \warning
%>  Currently, only the restart output files in ASCII format can be read via this method.<br>
%>  The binary restart files are not meant to be parsed via this method.<br>
%>  To request for ASCII restart output files in simulations,
%>  set the input simulation specification,<br>
%>  \code{.m}
%>      sampler.spec.outputRestartFileFormat = "ascii";
%>  \endcode
%>  where ``sampler`` can be an instance of any one of the ParaMonte
%>  sampler classes, such as [pm.sampling.Paradram](@ref Paradram).<br>
%>
%>  \warning
%>  Avoid using this routine for very large long simulations.<br>
%>  Reading the full restart file of a large-scale simulation problem
%>  can be extremely memory-intensive.<br>
%>
%>  \note
%>  See the documentation of the subclasses of [pm.sampling.Sampler](@ref Sampler)
%>  (e.g., [pm.sampling.Paradram](@ref Paradram)) for example usage in action.<br>
%>
%>  \example{readRestart}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>
%>      restartList = pm.sampling.readRestart([], "./out/test_run_");
%>      restartList = pm.sampling.readRestart(sampler, "./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>
%>      restartList = pm.sampling.readRestart(sampler, []);
%>      restartList = pm.sampling.readRestart(sampler);
%>
%>  \endcode
%>
%>  \final{readRestart}
%>
%>  \author
%>  \AmirShahmoradi, Saturday Nov 2 2024, 11:11 PM, Dallas, TX.<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function restartList = readRestart(sampler, pattern)

    if  nargin < 1
        sampler = [];
    end
    if  isempty(sampler)
        sampler = pm.sampling.Sampler();
    end

    if  nargin < 2
        if  0 < pm.array.len(sampler.spec.outputFileName)
            pattern = string(sampler.spec.outputFileName);
        else
            pattern = [];
        end
    end

    pathList = sampler.findfile("restart", pattern);
    restartList = cell(numel(pathList), 1);
    for ifile = numel(pathList) : -1 : 1

        if ~sampler.silent
            disp("processing file: """ + pathList(ifile) + """");
        end

        if  isequal(class(sampler), "pm.sampling.Sampler")
            restartList{ifile} = pm.sampling.FileContentsRestart(pathList(ifile), sampler.silent);
        elseif isequal(class(sampler), "pm.sampling.Paradram")
            restartList{ifile} = pm.sampling.FileContentsRestartDRAM(pathList(ifile), sampler.silent);
        else
            help("pm.sampling.readRestart")
            disp("class(sampler)");
            disp( class(sampler) );
            error   ( newline ...
                    + "Invalid input ``sampler`` object type." + newline ...
                    + "See the documentation displayed above for more " + newline ...
                    + "information on possible values for ``sampler`` argument." + newline ...
                    + newline ...
                    );
        end

    end

end