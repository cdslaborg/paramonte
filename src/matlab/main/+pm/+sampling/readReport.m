%>  \brief
%>  Return a list of objects of superclass [pm.sampling.FileContentsReport](@ref FileContentsReport)
%>  containing the contents of a (set of) ParaMonte simulation output report file(s) whose paths match
%>  the specified input ``pattern`` or the simulation specification ``sampler.spec.outputFileName``.<br>
%>
%>  \warning
%>  This function is to be only used for post-processing of the output report file(s) of an already finished simulation.<br>
%>  Although possible, this method is NOT meant to be called by all processes in MPI-parallel simulations.<br>
%>
%>  \param[in]  sampler :   The input object of superclass [pm.sampling.Sampler](@ref Sampler)
%>                          whose type and properties determine the type of the output object(s).<br>
%>                          (**optional**. If empty, it is set to [pm.sampling.Sampler()](@ref Sampler).)
%>  \param[in]  pattern :   The input scalar MATLAB string containing the pattern matching
%>                          the desired report file(s) whose contents is to be read.<br>
%>                          The specified ``pattern`` only needs to partially identify
%>                          the name of the simulation to which the report file belongs.<br>
%>                          For example, specifying ``"./mydir/mysim"`` as input will
%>                          lead to a search for file(s) beginning with "mysim" and
%>                          ending with ``"_report.txt"`` inside the directory ``"./mydir/"``.<br>
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
%>  ``reportList``      :   The output MATLAB cell array of objects of superclass [pm.sampling.FileContentsReport](@ref FileContentsReport),
%>                          each of which corresponds to the contents of a unique ParaMonte sampler report file.<br>
%>                          <ol>
%>                              <li>    If the input argument ``sampler`` is of type [pm.sampling.Sampler](@ref Sampler),
%>                                      then the output array elements are of type [pm.sampling.FileContentsReport](@ref FileContentsReport).<br>
%>                              <li>    If the input argument ``sampler`` is of type [pm.sampling.Paradram](@ref Paradram),
%>                                      then the output array elements are of type [pm.sampling.FileContentsReport](@ref FileContentsReport).<br>
%>                          </ol>
%>
%>  \interface{readReport}
%>  \code{.m}
%>
%>      reportList = pm.sampling.readReport();
%>      reportList = pm.sampling.readReport([]);
%>      reportList = pm.sampling.readReport([], []);
%>      reportList = pm.sampling.readReport([], pattern);
%>      reportList = pm.sampling.readReport(sampler, []);
%>      reportList = pm.sampling.readReport(sampler, pattern);
%>
%>  \endcode
%>
%>  \note
%>  See the documentation of the sampler subclasses
%>  (e.g., [pm.sampling.Paradram](@ref Paradram)) for example usage in action.<br>
%>
%>  \example{readReport}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>
%>      reportList = pm.sampling.readReport("./out/test_run_");
%>      reportList = pm.sampling.readReport([], "./out/test_run_");
%>      reportList = pm.sampling.readReport(sampler, "./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>
%>      reportList = pm.sampling.readReport(sampler, []);
%>      reportList = pm.sampling.readReport(sampler);
%>
%>  \endcode
%>
%>  \final{readReport}
%>
%>  \author
%>  \AmirShahmoradi, Saturday Nov 2 2024, 11:11 PM, Dallas, TX.<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function reportList = readReport(sampler, pattern)

    if  nargin < 2
        if  0 < pm.array.len(sampler.spec.outputFileName)
            pattern = string(sampler.spec.outputFileName);
        else
            pattern = [];
        end
    end

    if  nargin < 1
        sampler = [];
    end
    if  isempty(sampler)
        sampler = pm.sampling.Sampler();
    end

    pathList = sampler.findfile("report", pattern);
    reportList = cell(numel(pathList), 1);

    for ifile = numel(pathList) : -1 : 1

        if ~sampler.silent
            disp("processing file: """ + pathList(ifile) + """");
        end

        if  isequal(class(sampler), "pm.sampling.Sampler") || isequal(class(sampler), "pm.sampling.Paradram")
            reportList{ifile} = pm.sampling.FileContentsReport(pathList(ifile), sampler.silent, sampler.method);
        else
            help("pm.sampling.readReport")
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