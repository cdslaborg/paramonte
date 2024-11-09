%>  \brief
%>  Return a list of objects of superclass [pm.sampling.FileContentsChain](@ref FileContentsChain)
%>  containing the contents of a (set of) ParaMonte simulation output chain file(s) whose paths match
%>  the specified input ``pattern`` or the simulation specification ``sampler.spec.outputFileName``.<br>
%>
%>  \warning
%>  This function is to be only used for post-processing of the output chain file(s) of an already finished simulation.<br>
%>  Although possible, this method is NOT meant to be called by all processes in MPI-parallel simulations.<br>
%>
%>  \param[in]  sampler     :   The input object of superclass [pm.sampling.Sampler](@ref Sampler)
%>                              whose type and properties determine the type of the output object(s).<br>
%>                              (**optional**. If empty, it is set to [pm.sampling.Sampler()](@ref Sampler).)
%>  \param[in]  pattern     :   The input scalar MATLAB string containing the pattern matching
%>                              the desired chain file(s) whose contents is to be read.<br>
%>                              The specified ``pattern`` only needs to partially identify
%>                              the name of the simulation to which the chain file belongs.<br>
%>                              For example, specifying ``"./mydir/mysim"`` as input will
%>                              lead to a search for file(s) beginning with "mysim" and
%>                              ending with ``"_chain.txt"`` inside the directory ``"./mydir/"``.<br>
%>                              If there are multiple files matching in the input ``pattern``,
%>                              then all such files will be read and returned as elements of a list.<br>
%>                              If the specified pattern is a valid existing URL, the file will be
%>                              downloaded as a temporary file to the local system, its contents
%>                              shall be parsed and the file will be subsequently removed.<br>
%>                              If the input ``pattern`` is empty, then the method will search
%>                              for any possible candidate files with the appropriate suffix
%>                              in the current working directory.<br>
%>                              (**optional**. If empty, it is set to ``sampler.spec.outputFileName`` or if empty, it is set to ``"./"``.)
%>  \param[in]  sep         :   The input MATLAB string containing the field separator used in the file(s).<br>
%>                              (**optional**. If empty, it is set to ``sampler.spec.outputSeparator`` or if empty, it is automatically inferred.)
%>  \param[in]  varargin    :   Any set of arguments that must be directly passed to the relevant file contents parser class.<br>
%>                              (**optional**. default is empty, corresponding to no extra input arguments.)
%>
%>  \return
%>  ``chainList``           :   The output MATLAB cell array of objects of superclass [pm.sampling.FileContentsChain](@ref FileContentsChain),
%>                              each of which corresponds to the contents of a unique ParaMonte sampler chain file.<br>
%>                              <ol>
%>                                  <li>    If the input argument ``sampler`` is of type [pm.sampling.Sampler](@ref Sampler),
%>                                          then the output array elements are of type [pm.sampling.FileContentsChain](@ref FileContentsChain).<br>
%>                                  <li>    If the input argument ``sampler`` is of type [pm.sampling.Paradram](@ref Paradram),
%>                                          then the output array elements are of type [pm.sampling.FileContentsChainDRAM](@ref FileContentsChainDRAM).<br>
%>                              </ol>
%>
%>  \interface{readChain}
%>  \code{.m}
%>
%>      chainList = pm.sampling.readChain();
%>      chainList = pm.sampling.readChain([]);
%>      chainList = pm.sampling.readChain([], []);
%>      chainList = pm.sampling.readChain([], [], []);
%>      chainList = pm.sampling.readChain([], pattern);
%>      chainList = pm.sampling.readChain([], [], sep);
%>      chainList = pm.sampling.readChain([], pattern, sep);
%>      chainList = pm.sampling.readChain(sampler, pattern, sep, varargin);
%>
%>  \endcode
%>
%>  \note
%>  See the documentation of the ParaMonte MATLAB samplers
%>  (e.g., [pm.sampling.Paradram](@ref Paradram)) for example usage.<br>
%>
%>  \example{readChain}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>
%>      chainList = pm.sampling.readChain([], "./out/test_run_");
%>      chainList = pm.sampling.readChain(sampler, "./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>
%>      chainList = pm.sampling.readChain(sampler, [], ",");
%>      chainList = pm.sampling.readChain(sampler);
%>
%>      chainList = pm.sampling.readChain("./out/test_run_", ",");
%>
%>      sampler.spec.outputSeparator = ",";
%>
%>      chainList = pm.sampling.readChain(sampler, "./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.spec.outputSeparator = ",";
%>
%>      chainList = pm.sampling.readChain(sampler);
%>
%>  \endcode
%>
%>  \final{readChain}
%>
%>  \author
%>  \AmirShahmoradi, Saturday Nov 2 2024, 10:32 PM, Dallas, TX.<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
function chainList = readChain(sampler, pattern, sep, varargin)

    if  nargin < 3
        if  0 < pm.array.len(sampler.spec.outputSeparator)
            sep = string(sampler.spec.outputSeparator);
        else
            sep = [];
        end
    end

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

    pathList = sampler.findfile("chain", pattern);
    chainList = cell(numel(pathList), 1);

    for ifile = numel(pathList) : -1 : 1

        if ~sampler.silent
            disp("processing file: """ + pathList(ifile) + """");
        end

        if  isequal(class(sampler), "pm.sampling.Sampler")
            chainList{ifile} = pm.sampling.FileContentsChain(pathList(ifile), sampler.silent, sep, varargin{:});
        elseif isequal(class(sampler), "pm.sampling.Paradram")
            chainList{ifile} = pm.sampling.FileContentsChainDRAM(pathList(ifile), sampler.silent, sep, varargin{:});
        else
            help("pm.sampling.readChain")
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