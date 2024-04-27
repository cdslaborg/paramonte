function chainMarkovList = readChainMarkov(self, pattern, sep)
    %
    %   Return a list of objects of class ``pm.sampling.FileContentsChain``
    %   containing the content(s) of the ParaMonte simulation output chain
    %   file(s) whose path(s) match the specified input ``pattern`` or the
    %   simulation specification ``sampler.spec.outputFileName``.
    %
    %       \warning
    %
    %           Avoid using this routine for very large compact chains.
    %           Reading the full Markov chain of large-scale simulation problems
    %           can be extremely memory-intensive without any potential benefits.
    %
    %       \warning
    %
    %           This method is to be only used for post-processing of the output
    %           chain file(s) of an already finished simulation. It is NOT meant to
    %           be called by all processes in parallel mode, although it is possible.
    %
    %       \note
    %
    %           This routine is identical to ``readChain()`` method, except for the
    %           fact that upon reading the output chain files, it will also convert
    %           the chain contents from the default efficient compact format stored
    %           in the file to the full verbose Markov chain format.
    %
    %   Parameters
    %   ----------
    %
    %       pattern
    %
    %           See the documentation of the corresponding argument of
    %           the constructor of the method ``pm.sampling.Sampler.readChain``.
    %           (**optional**, The default is set by ``pm.sampling.Sampler.readChain``)
    %
    %       sep
    %
    %           See the documentation of the corresponding argument of
    %           the constructor of the method ``pm.sampling.Sampler.readChain``.
    %           (**optional**, The default is set by ``pm.sampling.Sampler.readChain``)
    %
    %   Returns
    %   -------
    %
    %       chainMarkovList
    %
    %           A cell array of objects of class ``pm.sampling.FileContentsChain``,
    %           each of which corresponds to the contents of a unique chain file.
    %           Try ``doc pm.sampling.FileContentsChain`` to see the documentation
    %           of the contents of the objects of the output list.
    %
    %   Interface
    %   ---------
    %
    %       sampler = pm.sampling.Paradram();
    %       sampler.readChainMarkov()
    %       sampler.readChainMarkov([])
    %       sampler.readChainMarkov(file)
    %       sampler.readChainMarkov([], [])
    %       sampler.readChainMarkov(file, [])
    %       sampler.readChainMarkov(file, sep)
    %
    %   Example
    %   -------
    %
    %       sampler = pm.sampling.Paradram()
    %
    %       sampler.readChainMarkov();
    %
    %       sim.readChainMarkov("./out/test_run_");
    %
    %       sim.spec.outputFileName = "./out/test_run_";
    %       sim.readChainMarkov();
    %
    %       pmpd.readChainMarkov("./out/test_run_", ",");
    %
    %       pmpd.spec.outputSeparator = ",";
    %       pmpd.readChainMarkov("./out/test_run_");
    %
    %       pmpd.spec.outputFileName = "./out/test_run_";
    %       pmpd.spec.outputSeparator = ",";
    %       pmpd.readChainMarkov();
    %
    %       Except the first ``readChainMarkov()`` example,
    %       all other examples are functionally equivalent.
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
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