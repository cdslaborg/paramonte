classdef SpecDRAM < pm.sampling.SpecMCMC
    %
    %   This is the base class for the ParaMonte sampler DRAM specifications.
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
    %       The DRAM simulation specifications are all described on this page:
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
    %       An object of class ``pm.sampling.SpecDRAM``.
    %
    %   Interface
    %   ---------
    %
    %       spec = pm.sampling.SpecDRAM()
    %       spec = pm.sampling.SpecDRAM([])
    %       spec = pm.sampling.SpecDRAM([], [])
    %       spec = pm.sampling.SpecDRAM(method)
    %       spec = pm.sampling.SpecDRAM(method, [])
    %       spec = pm.sampling.SpecDRAM([], silent)
    %       spec = pm.sampling.SpecDRAM(method, silent)
    %
    %   LICENSE
    %   -------
    %
    %       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
    %
    properties
        proposalDelayedRejectionScaleFactor = [];
        proposalDelayedRejectionCount       = [];
        proposalAdaptationPeriod            = [];
        proposalAdaptationCountGreedy       = [];
        proposalAdaptationCount             = [];
        burninAdaptationMeasure             = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        function self = SpecDRAM(method, silent)
            %
            %   Construct and return an object of class ``pm.sampling.SpecDRAM``.
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
            %       The output scalar object of class ``pm.sampling.SpecDRAM``.
            %
            %   Interface
            %   ---------
            %
            %       spec = pm.sampling.SpecDRAM()
            %       spec = pm.sampling.SpecDRAM([])
            %       spec = pm.sampling.SpecDRAM([], [])
            %       spec = pm.sampling.SpecDRAM(method)
            %       spec = pm.sampling.SpecDRAM(method, [])
            %       spec = pm.sampling.SpecDRAM([], silent)
            %       spec = pm.sampling.SpecDRAM(method, silent)
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
                method = "ParaDRAM";
            end
            self = self@pm.sampling.SpecMCMC(method, silent);
            self.url = "https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/";
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
            entries = getEntriesNML@pm.sampling.SpecMCMC(self, ndim);
            if  ~isempty(self.proposalDelayedRejectionScaleFactor); entries = entries + self.nmlsep + pm.introspection.getEntryNML("burninAdaptationMeasure              ", self.burninAdaptationMeasure            , "real"   , 1); end
            if  ~isempty(self.proposalDelayedRejectionCount      ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalAdaptationCount              ", self.proposalAdaptationCount            , "integer", 1); end
            if  ~isempty(self.proposalAdaptationPeriod           ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalAdaptationCountGreedy        ", self.proposalAdaptationCountGreedy      , "integer", 1); end
            if  ~isempty(self.proposalAdaptationCountGreedy      ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalAdaptationPeriod             ", self.proposalAdaptationPeriod           , "integer", 1); end
            if  ~isempty(self.proposalAdaptationCount            ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalDelayedRejectionCount        ", self.proposalDelayedRejectionCount      , "integer", 1); end
            if  ~isempty(self.burninAdaptationMeasure            ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalDelayedRejectionScaleFactor  ", self.proposalDelayedRejectionScaleFactor, "real"   , self.proposalDelayedRejectionCount); end
        end

    end

end