%>  \brief
%>  Return a list of objects of superclass [pm.sampling.FileContentsChain](@ref FileContentsChain)
%>  containing the content(s) of the ParaMonte simulation output chain file(s) whose path(s)
%>  match the specified input ``pattern`` or the simulation
%>  specification ``sampler.spec.outputFileName``.
%>
%>  \warning
%>  This method is to be only used for post-processing of the output chain file(s) of an already finished simulation.<br>
%>  Although possible, this method is NOT meant to be called by all processes in MPI-parallel simulations.<br>
%>
%>  \param[in]  self        :   The input parent object of superclass [pm.sampling.Sampler](@ref Sampler)
%>                              which is **implicitly** passed to this dynamic method (not by the user).<br>
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
%>                              (**optional**, default = ``sampler.spec.outputFileName`` or ``"./"``)
%>  \param[in]  sep         :   The input MATLAB string containing the field separator used in the file(s).<br>
%>                              (**optional**, default = ``sampler.spec.outputSeparator`` or automatically inferred.)
%>  \param[in]  varargin    :   Any set of arguments that must be directly passed to the relevant file contents parser class.<br>
%>                              (**optional**. default is empty, corresponding to no extra input arguments.)
%>
%>  \return
%>  ``chainList``           :   The output MATLAB cell array of objects
%>                              of superclass [pm.sampling.FileContentsChain](@ref FileContentsChain),
%>                              each of which corresponds to the contents of a unique chain file.<br>
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
%>      chainList = sampler.readChain(pattern, sep, varargin);
%>
%>  \endcode
%>
%>  \note
%>  See the documentation of the subclasses of [pm.sampling.Sampler](@ref Sampler)
%>  (e.g., [pm.sampling.Paradram](@ref Paradram)) for example usage in action.<br>
%>
%>  \example{readChain}
%>  \code{.m}
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
%>  \endcode
%>
%>  \final{readChain}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:18 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function chainList = readChain(self, pattern, sep, varargin)

    if  nargin < 3
        sep = [];
    end

    if  nargin < 2
        pattern = [];
    end

    chainList = pm.sampling.readChain(self, pattern, sep, varargin{:});

end