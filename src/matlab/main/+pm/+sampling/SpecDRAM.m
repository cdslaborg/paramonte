%>  \brief
%>  This is the base class for the ParaMonte sampler DRAM specifications.<br>
%>
%>  \details
%>  This is an abstract class that is not meant to be used by the user.<br>
%>  See the class constructor documentation.<br>
%>
%>  \note
%>  All class attributes can be set after constructing an instance of this class.<br>
%>
%>  \note
%>  The DRAM simulation specifications are all described on this page:
%>  [ParaDRAM simulation specifications listing](\pmdoc_usage_sampling/paradram/specifications/)<br>
%>
%>  \note
%>  See below for information on the methods.<br>
%>
%>  \final{SpecDRAM}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 3:45 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef SpecDRAM < pm.sampling.SpecMCMC

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        %> \specdram{proposaladaptationburnin}
        proposalAdaptationBurnin = [];
        %> \specdram{proposaladaptationcount}
        proposalAdaptationCount = [];
        %> \specdram{proposaladaptationcountgreedy}
        proposalAdaptationCountGreedy = [];
        %> \specdram{proposaladaptationperiod}
        proposalAdaptationPeriod = [];
        %> \specdram{proposaldelayedrejectioncount}
        proposalDelayedRejectionCount = [];
        %> \specdram{proposaldelayedrejectionscale}
        proposalDelayedRejectionScale = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %>  \brief
        %>  Construct and return an object of class [pm.sampling.SpecDRAM](@ref SpecDRAM).<br>
        %>
        %>  \param[in]  method  :   The input scalar MATLAB string containing
        %>                          the name of the specific ParaMonte sampler
        %>                          whose simulation specifications are to be
        %>                          stored in the output of this constructor.<br>
        %>  \param[in]  silent  :   The input scalar MATLAB logical.<br>
        %>                          If ``true``, all descriptive messages on
        %>                          the MATLAB command line will be suppressed.<br>
        %>                          (**optional**, default = ``false``)
        %>
        %>  \return
        %>  The output scalar object of class [pm.sampling.SpecDRAM](@ref SpecDRAM).<br>
        %>
        %>  \interface{SpecDRAM}
        %>  \code{.m}
        %>
        %>      spec = pm.sampling.SpecDRAM()
        %>      spec = pm.sampling.SpecDRAM([])
        %>      spec = pm.sampling.SpecDRAM([], [])
        %>      spec = pm.sampling.SpecDRAM(method)
        %>      spec = pm.sampling.SpecDRAM(method, [])
        %>      spec = pm.sampling.SpecDRAM([], silent)
        %>      spec = pm.sampling.SpecDRAM(method, silent)
        %>
        %>  \endcode
        %>
        %>  \final{SpecDRAM}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 3:46 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SpecDRAM(method, silent)
            if  nargin < 2
                silent = [];
            end
            if  nargin < 1
                method = "ParaDRAM";
            end
            self = self@pm.sampling.SpecMCMC(method, silent);
            self.url = "https://www.cdslab.org/paramonte/generic/" + pm.lib.version("generic", "major") + "/usage/sampling/paradram/specifications/";
        end

    end

    methods(Hidden)

        %>  \brief
        %>  Ensure all specification properties of the parent object are sensible.<br>
        %>
        %>  \details
        %>  This is a dynamic method of the class [pm.sampling.SpecDRAM](@ref SpecDRAM).<br>
        %>
        %>  \param[in]  self    :   The input parent object of class [pm.sampling.SpecDRAM](@ref SpecDRAM)
        %>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
        %>  \param[in]  ndim    :   The input scalar MATLAB integer containing
        %>                          the number of dimensions of the domain of the
        %>                          object function that is to be explored.<br>
        %>
        %>  \return
        %>  ``entries``         :   The output scalar MATLAB string containing
        %>                          the simulation specifications converted to
        %>                          a Fortran-namelist-compatible entry.<br>
        %>
        %>  \interface{getEntriesNML}
        %>  \code{.m}
        %>
        %>      entries = self.getEntriesNML(ndim)
        %>
        %>  \endcode
        %>
        %>  \final{getEntriesNML}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 3:48 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function entries = getEntriesNML(self, ndim)
            entries = getEntriesNML@pm.sampling.SpecMCMC(self, ndim);
            if  ~isempty(self.proposalAdaptationBurnin      ); entries = entries + self.nmlsep + pm.fortran.getEntryNML("proposalAdaptationBurnin     ", self.proposalAdaptationBurnin        , "real"   , 1); end
            if  ~isempty(self.proposalAdaptationCount       ); entries = entries + self.nmlsep + pm.fortran.getEntryNML("proposalAdaptationCount      ", self.proposalAdaptationCount         , "integer", 1); end
            if  ~isempty(self.proposalAdaptationCountGreedy ); entries = entries + self.nmlsep + pm.fortran.getEntryNML("proposalAdaptationCountGreedy", self.proposalAdaptationCountGreedy   , "integer", 1); end
            if  ~isempty(self.proposalAdaptationPeriod      ); entries = entries + self.nmlsep + pm.fortran.getEntryNML("proposalAdaptationPeriod     ", self.proposalAdaptationPeriod        , "integer", 1); end
            if  ~isempty(self.proposalDelayedRejectionCount ); entries = entries + self.nmlsep + pm.fortran.getEntryNML("proposalDelayedRejectionCount", self.proposalDelayedRejectionCount   , "integer", 1); end
            if  ~isempty(self.proposalDelayedRejectionScale ); entries = entries + self.nmlsep + pm.fortran.getEntryNML("proposalDelayedRejectionScale", self.proposalDelayedRejectionScale   , "real"   , intmax); end
        end

    end

end