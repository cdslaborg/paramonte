classdef SpecDRAM_ProposalModel_class < handle

    properties (Constant)
        CLASS_NAME              = "@SpecDRAM_ProposalModel_class"
        MAX_LEN_PROPOSAL_MODEL  = 63
    end

    properties
        isUniform               = []
        isNormal                = []
        uniform                 = []
        normal                  = []
        val                     = []
        def                     = []
        desc                    = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecDRAM_ProposalModel_class()
            self.isUniform  = false;
            self.isNormal   = false;
            self.uniform    = "uniform";
            self.normal     = "normal";
            self.def        = self.normal;
            self.desc       = "proposalModel is a string variable containing the name of the proposal distribution for the MCMC sampler. "              ...
                            + "The string value must be enclosed by either single or double quotation marks when provided as input. "                   ...
                            + "One option is currently supported:" + newline + newline                                                                  ...
                            + "    proposalModel = '" + self.uniform + "'" + newline + newline                                                          ...
                            + "            The proposals will be drawn uniformly from within a ndim-dimensional ellipsoid whose covariance matrix "     ...
                            + "and scale are initialized by the user and optionally adaptively updated throughout the simulation." + newline + newline  ...
                            + "The default value is '" + self.def + "'."                                                                                ...
                            ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, proposalModel)
            if isempty(proposalModel)
                self.val = self.def;
            else
                self.val = lower(strtrim(proposalModel));
            end
            if self.val == self.normal,     self.isNormal   = true; end
            if self.val == self.uniform,    self.isUniform  = true; end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if ~(self.isNormal || self.isUniform)
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                        ...
                                + "Invalid requested value for the proposalModel of " + methodName + ". The input requested "   ...
                                + "proposal model (" + self.val + ") is not supported. "                                        ...
                                + "The variable proposalModel cannot be set to anything other than '"                           ...
                                + self.normal + "', or '" + self.uniform + "'." + newline + newline                             ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end