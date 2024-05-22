%>  \brief
%>  Return a list of objects of class ``pm.sampling.FileContentsChain``
%>  containing the content(s) of the ParaMonte simulation output chain
%>  file(s) whose path(s) match the specified input ``pattern`` or the
%>  simulation specification ``sampler.spec.outputFileName``.
%>
%>  \warning
%>  This method is to be only used for post-processing of the output
%>  chain file(s) of an already finished simulation. Although possible,
%>  this method is NOT meant to be called by all processes
%>  in MPI-parallel simulations.
%>
%>  \param[in]  pattern :   The input scalar MATLAB string containing the pattern matching
%>                          the desired chain file(s) whose contents is to be read.
%>                          The specified ``pattern`` only needs to partially identify
%>                          the name of the simulation to which the chain file belongs.
%>                          For example, specifying ``"./mydir/mysim"`` as input will
%>                          lead to a search for file(s) beginning with "mysim" and
%>                          ending with ``"_chain.txt"`` inside the directory ``"./mydir/"``.<br>                    
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
%>  \param[in]  sep     :   The input MATLAB string containing the field separator used in the file(s).
%>                          (optional, default = ``sampler.spec.outputSeparator`` or automatically inferred.)
%>
%>  \return
%>  `chainList`         :   The output MATLAB cell array of objects
%>                          of class ``pm.sampling.FileContentsChain``,
%>                          each of which corresponds to the contents
%>                          of a unique chain file.
%>
%>  \interface{readChain}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>      chainList = sampler.readChain();
%>      chainList = sampler.readChain([]);
%>      chainList = sampler.readChain([], []);
%>      chainList = sampler.readChain(pattern);
%>      chainList = sampler.readChain([], sep);
%>      chainList = sampler.readChain(pattern, sep);
%>
%>  \endcode
%>
%>  \example{readChain}
%>
%>      sampler.readChain("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.readChain();
%>
%>      sampler.readChain("./out/test_run_", ",");
%>
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readChain("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readChain();
%>
%>  \final{readChain}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:18 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
%>
function chainList = readChain(self, pattern, sep)
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
    ftype = "chain";
    pathList = self.findfile(ftype, pattern);
    chainList = cell(length(pathList), 1);
    for ifile = length(pathList) : -1 : 1
        if ~self.silent
            disp("processing file: """ + pathList(ifile) + """");
        end
        chainList{ifile} = pm.sampling.FileContentsChain(pathList(ifile), self.silent, sep);
    end
end