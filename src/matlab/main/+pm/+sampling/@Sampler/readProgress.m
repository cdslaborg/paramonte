%>  \brief
%>  Return a list of objects of class [pm.sampling.FileContentsProgress](@ref FileContentsProgress)
%>  containing the content(s) of the ParaMonte simulation output progress
%>  file(s) whose path(s) match the specified input ``pattern`` or the
%>  simulation specification ``sampler.spec.outputFileName``.<br>
%>
%>  \warning
%>  This method is to be only used for post-processing of the output
%>  progress file(s) of an already finished simulation. Although possible,
%>  this method is NOT meant to be called by all processes
%>  in MPI-parallel simulations.<br>
%>
%>  \param[in]  pattern :   The input scalar MATLAB string containing the pattern matching
%>                          the desired progress file(s) whose contents is to be read.<br>
%>                          The specified ``pattern`` only needs to partially identify
%>                          the name of the simulation to which the progress file belongs.<br>
%>                          For example, specifying ``"./mydir/mysim"`` as input will
%>                          lead to a search for file(s) beginning with "mysim" and
%>                          ending with ``"_progress.txt"`` inside the directory ``"./mydir/"``.<br>
%>                          If there are multiple files matching in the input ``pattern``,
%>                          then all such files will be read and returned as elements of a list.<br>
%>                          If the specified pattern is a valid existing URL, the file will be
%>                          downloaded as a temporary file to the local system, its contents
%>                          shall be parsed and the file will be subsequently removed.<br>
%>                          If the input ``pattern`` is empty, then the method will search
%>                          for any possible candidate files with the appropriate suffix
%>                          in the current working directory.<br>
%>                          (optional, default = ``sampler.spec.outputFileName`` or ``"./"``)
%>  \param[in]  sep     :   The input MATLAB string containing the field separator used in the file(s).<br>
%>                          (optional, default = ``sampler.spec.outputSeparator`` or automatically inferred.)
%>
%>  \return
%>  ``progressList``    :   The output MATLAB cell array of objects
%>                          of class [pm.sampling.FileContentsProgress](@ref FileContentsProgress),
%>                          each of which corresponds to the contents
%>                          of a unique progress file.<br>
%>
%>  \interface{readProgress}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>      progressList = sampler.readProgress();
%>      progressList = sampler.readProgress([]);
%>      progressList = sampler.readProgress([], []);
%>      progressList = sampler.readProgress(pattern);
%>      progressList = sampler.readProgress([], sep);
%>      progressList = sampler.readProgress(pattern, sep);
%>
%>  \endcode
%>
%>  \final{readProgress}
%>
%>  \note
%>  See the documentation of the sampler subclasses
%>  (e.g., [pm.sampling.Paradram](@ref Paradram)) for example usage in action.<br>
%>
%>  \example{readProgress}
%>  \code{.m}
%>
%>      sampler.readProgress("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.readProgress();
%>
%>      sampler.readProgress("./out/test_run_", ",");
%>
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readProgress("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readProgress();
%>
%>  \endcode
%>
%>  \final{readProgress}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:25 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function progressList = readProgress(self, pattern, sep)
    if nargin < 3
        if 0 < pm.array.len(self.spec.outputSeparator)
            sep = string(self.spec.outputSeparator);
        else
            sep = [];
        end
    end
    if nargin < 2
        if 0 < pm.array.len(self.spec.outputFileName)
            pattern = string(self.spec.outputFileName);
        else
            pattern = [];
        end
    end
    ftype = "progress";
    pathList = self.findfile(ftype, pattern);
    progressList = cell(length(pathList), 1);
    for ifile = length(pathList) : -1 : 1
        if ~self.silent
            disp("processing file: """ + pathList(ifile) + """");
        end
        progressList{ifile} = pm.sampling.FileContentsProgress(pathList(ifile), self.silent, sep);
    end
end