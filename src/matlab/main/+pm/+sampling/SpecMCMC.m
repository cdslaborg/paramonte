%>  \brief
%>  This is the base class for the ParaMonte sampler MCMC specifications.
%>  This is an abstract class that is not meant to be used by the user.
%>  See the class constructor.
%>
%>  \note
%>  All class attributes can be set after constructing an instance of this class.
%>
%>  \note
%>  The MCMC simulation specifications are all described on this page:
%>  https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/
%>
%>  \note
%>  See below for information on the methods.
%>
%>  \return
%>  An object of class ``pm.sampling.SpecMCMC``.
%>
%>  \interface{SpecMCMC}
%>  \code{.m}
%>
%>      spec = pm.sampling.SpecMCMC()
%>      spec = pm.sampling.SpecMCMC([])
%>      spec = pm.sampling.SpecMCMC([], [])
%>      spec = pm.sampling.SpecMCMC(method)
%>      spec = pm.sampling.SpecMCMC(method, [])
%>      spec = pm.sampling.SpecMCMC([], silent)
%>      spec = pm.sampling.SpecMCMC(method, silent)
%>
%>  \endcode
%>  \final{SpecMCMC}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:54 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
classdef SpecMCMC < pm.sampling.SpecBase

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        outputChainSize                     = [];
        outputSampleRefinementCount         = [];
        outputSampleRefinementMethod        = [];
        proposal                            = [];
        proposalCor                         = [];
        proposalCov                         = [];
        proposalScale                       = [];
        proposalStart                       = [];
        proposalStartDomainCubeLimitLower   = [];
        proposalStartDomainCubeLimitUpper   = [];
        proposalStartRandomized             = [];
        proposalStd                         = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)
        %>  \brief
        %>  Construct and return an object of class ``pm.sampling.SpecMCMC``.
        %>
        %>  \param[in]  method  :   The input scalar MATLAB string containing
        %>                          the name of the specific ParaMonte sampler
        %>                          whose simulation specifications are to be
        %>                          stored in the output of this constructor.
        %>  
        %>  \param[in]  silent  :   The input scalar MATLAB logical.
        %>                          If ``true``, all descriptive messages on
        %>                          the MATLAB command line will be suppressed.
        %>                          (**optional**, default = ``false``)
        %>  
        %>  \return
        %>  The output scalar object of class ``pm.sampling.SpecMCMC``.
        %>
        %>  \interface{SpecMCMC}
        %>  \code{.m}
        %>
        %>      spec = pm.sampling.SpecMCMC()
        %>      spec = pm.sampling.SpecMCMC([])
        %>      spec = pm.sampling.SpecMCMC([], [])
        %>      spec = pm.sampling.SpecMCMC(method)
        %>      spec = pm.sampling.SpecMCMC(method, [])
        %>      spec = pm.sampling.SpecMCMC([], silent)
        %>      spec = pm.sampling.SpecMCMC(method, silent)
        %>
        %>  \endcode
        %>  \final{SpecMCMC}
        %>
        %>  \author
        %>  \JoshuaOsborne, May 21 2024, 12:58 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function self = SpecMCMC(method, silent)
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
        %>  \brief
        %>  Ensure all specification properties of the parent object are sensible.
        %>  This is a dynamic method of the class ``pm.sampling.SpecMCMC``.
        %>
        %>  \param[in]  ndim    :   The input scalar MATLAB integer containing the
        %>                          number of dimensions of the domain of the
        %>                          object function that is to be explored.
        %>
        %>  \return
        %>  `entries`           :   The output scalar MATLAB string containing
        %>                          the simulation specifications converted to
        %>                          a Fortran-namelist-compatible entry.
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
        %>  \JoshuaOsborne, May 21 2024, 1:01 AM, University of Texas at Arlington<br>
        %>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center, Washington, D.C.<br>
        %>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
        function entries = getEntriesNML(self, ndim)
            entries = getEntriesNML@pm.sampling.SpecBase(self, ndim);
            if  ~isempty(self.outputChainSize                   ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("outputChainSize                  ", self.outputChainSize                     , "integer" , 1); end
            if  ~isempty(self.outputSampleRefinementCount       ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("outputSampleRefinementCount      ", self.outputSampleRefinementCount         , "integer" , 1); end
            if  ~isempty(self.outputSampleRefinementMethod      ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("outputSampleRefinementMethod     ", self.outputSampleRefinementMethod        , "string"  , 1); end
            if  ~isempty(self.proposal                          ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposal                         ", self.proposal                            , "string"  , 1); end
            if  ~isempty(self.proposalCor                       ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalCor                      ", self.proposalCor                         , "real"    , ndim^2); end
            if  ~isempty(self.proposalCov                       ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalCov                      ", self.proposalCov                         , "real"    , ndim^2); end
            if  ~isempty(self.proposalScale                     ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalScale                    ", self.proposalScale                       , "string"  , 1); end
            if  ~isempty(self.proposalStart                     ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStart                    ", self.proposalStart                       , "real"    , ndim); end
            if  ~isempty(self.proposalStartDomainCubeLimitLower ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStartDomainCubeLimitLower", self.proposalStartDomainCubeLimitLower   , "real"    , ndim); end
            if  ~isempty(self.proposalStartDomainCubeLimitUpper ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStartDomainCubeLimitUpper", self.proposalStartDomainCubeLimitUpper   , "real"    , ndim); end
            if  ~isempty(self.proposalStartRandomized           ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStartRandomized          ", self.proposalStartRandomized             , "logical" , 1); end
            if  ~isempty(self.proposalStd                       ); entries = entries + self.nmlsep + pm.introspection.getEntryNML("proposalStd                      ", self.proposalStd                         , "real"    , ndim); end
        end

    end

end