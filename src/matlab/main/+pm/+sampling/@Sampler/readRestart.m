%>  \brief
%>  Return a list of objects of class [pm.sampling.FileContentsRestartDRAM](@ref FileContentsRestartDRAM)
%>  containing the content(s) of the ParaMonte simulation output restart
%>  file(s) whose path(s) match the specified input ``pattern`` or the
%>  simulation specification ``sampler.spec.outputFileName``.<br>
%>
%>  \warning
%>  This method is to be only used for post-processing of the output
%>  restart file(s) of an already finished simulation.<br>
%>  Although possible, this method is NOT meant
%>  to be called by all processes in
%>  MPI-parallel simulations.<br>
%>
%>  \warning
%>  Currently, only the restart output files in ASCII format can be read via this method.<br>
%>  The binary restart files are not meant to be parsed via this method.<br>
%>  To request for ASCII restart output files in simulations,
%>  set the input simulation specification,
%>
%>  \code{.m}
%>      sampler.spec.restartFileFormat = "ascii"
%>  \endcode
%>
%>  where ``sampler`` can be an instance of any one of the ParaMonte
%>  sampler classes, such as [pm.sampling.Paradram](@ref Paradram).<br>
%>
%>  \warning
%>  Avoid using this routine for very large long simulations.<br>
%>  Reading the full restart file of a large-scale simulation problem
%>  can be extremely memory-intensive.<br>
%>
%>  \param[in]  pattern :   The input scalar MATLAB string containing the pattern matching
%>                          the desired restart file(s) whose contents is to be read.<br>
%>                          The specified ``pattern`` only needs to partially identify
%>                          the name of the simulation to which the restart file belongs.<br>
%>                          For example, specifying ``"./mydir/mysim"`` as input will
%>                          lead to a search for file(s) beginning with "mysim" and
%>                          ending with ``"_restart.txt"`` inside the directory ``"./mydir/"``.<br>
%>                          If there are multiple files matching in the input ``pattern``,
%>                          then all such files will be read and returned as elements of a list.<br>
%>                          If the specified pattern is a valid existing URL, the file will be
%>                          downloaded as a temporary file to the local system, its contents
%>                          shall be parsed and the file will be subsequently removed.<br>
%>                          If the input ``pattern`` is empty, then the method will search
%>                          for any possible candidate files with the appropriate suffix
%>                          in the current working directory.<br>
%>                          (optional, default = ``sampler.spec.outputFileName`` or ``"./"``)
%>
%>  \return
%>  ``restartList``     :   The output MATLAB cell array of objects
%>                          of class [pm.sampling.FileContentsRestart](@ref FileContentsRestart),
%>                          each of which corresponds to the contents
%>                          of a unique restart file.<br>
%>
%>  \interface{readRestart}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>      restartList = sampler.readRestart();
%>      restartList = sampler.readRestart([]);
%>      restartList = sampler.readRestart(pattern);
%>
%>  \endcode
%>
%>  \note
%>  See the documentation of the sampler subclasses
%>  (e.g., [pm.sampling.Paradram](@ref Paradram)) for example usage in action.<br>
%>
%>  \example{readRestart}
%>  \code{.m}
%>
%>      sampler.readRestart("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.readRestart();
%>
%>      sampler.readRestart("./out/test_run_", ",");
%>
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readRestart("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readRestart();
%>
%>  \endcode
%>
%>  \final{readRestart}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:35 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
%>
function restartList = readRestart(self, pattern)
    if nargin < 2
        if 0 < pm.array.len(self.spec.outputFileName)
            pattern = string(self.spec.outputFileName);
        else
            pattern = [];
        end
    end
    ftype = "restart";
    pathList = self.findfile(ftype, pattern);
    restartList = cell(length(pathList), 1);
    for ifile = length(pathList) : -1 : 1
        if ~self.silent
            disp("processing file: """ + pathList(ifile) + """");
        end
        restartList{ifile} = pm.sampling.FileContentsRestartDRAM(pathList(ifile), self.silent);
    end
end