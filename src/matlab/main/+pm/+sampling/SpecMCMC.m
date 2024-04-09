classdef SpecMCMC < pm.sampling.SpecBase
    %
    %   This is the base class for the ParaMonte sampler MCMC specifications.
    %
    %   This is an abstract class that is not meant to be used by the user.
    %
    %
    %   Parameters
    %   ----------
    %
    %       See the class constructor.
    %
    %       \note
    %
    %           All class attributes can be set after constructing an instance of this class.
    %
    %   Attributes
    %   ----------
    %
    %       The MCMC simulation specifications are all described on this page:
    %
    %           https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/
    %
    %   Methods
    %   -------
    %
    %       See below for information on the methods.
    %
    %   Returns
    %   -------
    %
    %       An object of class ``pm.sampling.SpecMCMC``.
    %
    %   Interface
    %   ---------
    %
    %       spec = pm.sampling.SpecMCMC()
    %       spec = pm.sampling.SpecMCMC([])
    %       spec = pm.sampling.SpecMCMC([], [])
    %       spec = pm.sampling.SpecMCMC(method)
    %       spec = pm.sampling.SpecMCMC(method, [])
    %       spec = pm.sampling.SpecMCMC([], silent)
    %       spec = pm.sampling.SpecMCMC(method, silent)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties
        outputChainSize                     = [];
        outputSampleRefinementCount         = [];
        outputSampleRefinementMethod        = [];
        proposal                            = [];
        proposalCorMat                      = [];
        proposalCovMat                      = [];
        proposalScaleFactor                 = [];
        proposalStart                       = [];
        proposalStartDomainCubeLimitLower   = [];
        proposalStartDomainCubeLimitUpper   = [];
        proposalStartRandomized             = [];
        proposalStdVec                      = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        function self = SpecMCMC(method, silent)
            %
            %   Construct and return an object of class ``pm.sampling.SpecMCMC``.
            %
            %   Parameters
            %   ----------
            %
            %       method
            %
            %           The input scalar MATLAB string containing
            %           the name of the specific ParaMonte sampler
            %           whose simulation specifications are to be
            %           stored in the output of this constructor.
            %
            %       silent
            %
            %           The input scalar MATLAB logical.
            %           If ``true``, all descriptive messages on
            %           the MATLAB command line will be suppressed.
            %           (**optional**, default = ``false``)
            %
            %   Returns
            %   -------
            %
            %       The output scalar object of class ``pm.sampling.SpecMCMC``.
            %
            %   Interface
            %   ---------
            %
            %       spec = pm.sampling.SpecMCMC()
            %       spec = pm.sampling.SpecMCMC([])
            %       spec = pm.sampling.SpecMCMC([], [])
            %       spec = pm.sampling.SpecMCMC(method)
            %       spec = pm.sampling.SpecMCMC(method, [])
            %       spec = pm.sampling.SpecMCMC([], silent)
            %       spec = pm.sampling.SpecMCMC(method, silent)
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            if  nargin < 2
                silent = [];
            end
            if  nargin < 1
                method = "ParaMCMC";
            end
            self = self@pm.sampling.SpecBase(method, silent);
        end

    end

    methods(Hidden)

        function entries = getEntriesNML(self, ndim)
            %
            %   Ensure all specification properties of the parent object are sensible.
            %
            %   This is a dynamic method of the class ``pm.sampling.SpecMCMC``.
            %
            %   Parameters
            %   ----------
            %
            %       ndim
            %
            %           The input scalar MATLAB integer containing the
            %           number of dimensions of the domain of the
            %           object function that is to be explored.
            %
            %   Returns
            %   -------
            %
            %       entries
            %
            %           The output scalar MATLAB string containing
            %           the simulation specifications converted to
            %           a Fortran-namelist-compatible entry.
            %
            %   Interface
            %   ---------
            %
            %       entries = self.getEntriesNML(ndim)
            %
            %   LICENSE
            %   -------
            %
            %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
            %
            entries = getEntriesNML@pm.sampling.SpecBase(self, ndim);
            if  ~isempty(self.outputChainSize                   ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("outputChainSize                  ", self.outputChainSize                     , "integer", 1); end
            if  ~isempty(self.outputSampleRefinementCount       ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("outputSampleRefinementCount      ", self.outputSampleRefinementCount         , "integer", 1); end
            if  ~isempty(self.outputSampleRefinementMethod      ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("outputSampleRefinementMethod     ", self.outputSampleRefinementMethod        , "string" , 1); end
            if  ~isempty(self.proposal                          ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposal                         ", self.proposal                            , "string" , 1); end
            if  ~isempty(self.proposalCorMat                    ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalCorMat                   ", self.proposalCorMat                      , "real"   , ndim^2); end
            if  ~isempty(self.proposalCovMat                    ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalCovMat                   ", self.proposalCovMat                      , "real"   , ndim^2); end
            if  ~isempty(self.proposalScaleFactor               ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalScaleFactor              ", self.proposalScaleFactor                 , "string"  , 1); end
            if  ~isempty(self.proposalStart                     ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStart                    ", self.proposalStart                       , "real"   , ndim); end
            if  ~isempty(self.proposalStartDomainCubeLimitLower ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStartDomainCubeLimitLower", self.proposalStartDomainCubeLimitLower   , "real"   , ndim); end
            if  ~isempty(self.proposalStartDomainCubeLimitUpper ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStartDomainCubeLimitUpper", self.proposalStartDomainCubeLimitUpper   , "real"   , ndim); end
            if  ~isempty(self.proposalStartRandomized           ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStartRandomized          ", self.proposalStartRandomized             , "logical", 1); end
            if  ~isempty(self.proposalStdVec                    ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStdVec                   ", self.proposalStdVec                      , "real"   , ndim); end
        end

    end

end