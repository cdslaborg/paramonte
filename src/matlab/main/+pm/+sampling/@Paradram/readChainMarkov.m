%>  \brief
%>  Return a list of objects of class [pm.sampling.FileContentsChain](@ref FileContentsChain)
%>  containing the content(s) of the ParaMonte simulation output chain file(s) whose path(s)
%>  match the specified input ``pattern`` or the simulation specification 
%>  [sampler.spec.outputFileName](\pmdoc_usage_sampling/paradram/specifications/#outputfilename).
%>
%>  \param[in]  self    :   The input/output parent object of class [pm.sampling.Paradram](@ref Paradram)
%>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
%>  \param[in]  pattern :   See the documentation of the corresponding argument of
%>                          the constructor of the method [pm.sampling.Sampler.readChain](@ref Sampler::readChain).<br>
%>                          (**optional**, The default is set by [pm.sampling.Sampler.readChain](@ref Sampler::readChain))
%>  \param[in]  sep     :   See the documentation of the corresponding argument of
%>                          the constructor of the method [pm.sampling.Sampler.readChain](@ref Sampler::readChain).<br>
%>                          (**optional**, The default is set by [pm.sampling.Sampler.readChain](@ref Sampler::readChain))
%>
%>  \return
%>  ``chainMarkovList`` :   The output cell array of objects of class
%>                          [pm.sampling.FileContentsChain](@ref FileContentsChain),
%>                          each of which corresponds to the contents of a unique chain file.<br>
%>                          Try ``doc pm.sampling.FileContentsChain`` to see the documentation
%>                          of the contents of the objects of the output list.<br>
%>
%>  \interface{Paradram}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Paradram();
%>      chainMarkovList = sampler.readChainMarkov();
%>      chainMarkovList = sampler.readChainMarkov([]);
%>      chainMarkovList = sampler.readChainMarkov(file);
%>      chainMarkovList = sampler.readChainMarkov([], []);
%>      chainMarkovList = sampler.readChainMarkov(file, []);
%>      chainMarkovList = sampler.readChainMarkov(file, sep);
%>
%>  \endcode
%>
%>  \warning
%>  Avoid using this routine for very large compact chains.<br>
%>  Reading the full Markov chain of large-scale simulation problems
%>  can be extremely memory-intensive without any potential benefits.<br>
%>
%>  \warning
%>  This method is to be only used for post-processing of the output
%>  chain file(s) of an already finished simulation. It is NOT meant to
%>  be called by all processes in parallel mode, although it is possible.<br>
%>
%>  \note
%>  This routine is identical to [pm.sampling.Sampler.readChain](@ref Sampler::readChain)
%>   method, except for the fact that upon reading the output chain files, it will also
%>  convert the chain contents from the default efficient compact format stored
%>  in the file to the full verbose Markov chain format.<br>
%>
%>  \example{readChainMarkov}
%>  \include{lineno} example/sampling/Paradram/readChainMarkov/main.m
%>  \include{lineno} example/sampling/Paradram/readChainMarkov/main.out.m
%>  \vis{readChainMarkov}
%>  \image html example/sampling/Paradram/readChainMarkov/readChainMarkov.domain.png width=700
%>  <br><br>
%>  \image html example/sampling/Paradram/readChainMarkov/readChainMarkov.traceplot.png width=700
%>  <br><br>
%>  \image html example/sampling/Paradram/readChainMarkov/readChainMarkov.proposalAdaptation.png width=700
%>
%>  \final{readChainMarkov}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:06 AM, University of Texas at Arlington<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function chainMarkovList = readChainMarkov(self, pattern, sep)
    if nargin < 3
        sep = [];
    end
    if nargin < 2
        pattern = [];
    end
    chainMarkovList = self.readChain(pattern, sep);
    for ichain = 1 : length(chainMarkovList)
        chainMarkovList{ichain}.df = pm.array.verbose(chainMarkovList{ichain}.df, 1, chainMarkovList{ichain}.df.sampleWeight);
        chainMarkovList{ichain}.df.sampleWeight(:) = 1;
    end
end