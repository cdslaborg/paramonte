classdef SpecMCMC_SampleRefinementMethod_class < handle

    properties (Constant)
        CLASS_NAME                          = "@SpecMCMC_SampleRefinementMethod_class"
        BATCH_MEANS                         = "BatchMeans"
        MAX_CUMSUM_AUTOCORR                 = "MaxCumSumAutoCorr"
        LEN_MAX_CUMSUM_AUTOCORR             = length(SpecMCMC_SampleRefinementMethod_class.MAX_CUMSUM_AUTOCORR)
        LEN_BATCH_MEANS                     = length(SpecMCMC_SampleRefinementMethod_class.BATCH_MEANS)
        MAX_LEN_SAMPLE_REFINEMENT_METHOD    = 63
    end

    properties
        def                                 = []
        val                                 = []
        desc                                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecMCMC_SampleRefinementMethod_class(methodName)
            % self.isMaxCumSumAutoCorr    = false;
            % self.isViaCompactChain      = false;
            % self.isBatchMeans           = false;
            % self.isBatchMeans           = false;
            % self.isMaxCumSumAutoCorr    = false;
            % self.batchMeans             = "BatchMeans";
            % self.maxCumSumAutoCorr      = "MaxCumSumAutoCorr";
            % self.def                    = self.batchMeans;

            self.def                    = self.BATCH_MEANS;

            self.desc                   = "sampleRefinementMethod is a string variable that represents the method of computing the Integrated Autocorrelation Time "    ...
                                        + "(IAC) to be used in " + methodName + " for refining the final output MCMC chain and sample. "                                ...
                                        + "The string value must be enclosed by either single or double quotation marks when provided as input. "                       ...
                                        + "Options that are currently "                                                                                                 ...
                                        + "supported include:+" + newline + newline                                                                                     ...
                                        + "    sampleRefinementMethod = '" + self.BATCH_MEANS + "'" + newline + newline                                                 ...
                                        + "            This method of computing the Integrated Autocorrelation Time is based on the approach described in "             ...
                                                    + "SCHMEISER, B., 1982, Batch size effects in the analysis of simulation output, Oper. Res. 30 556-568. The "       ...
                                                    + "batch sizes in the BatchMeans method are chosen to be int(N^(2/3)) where N is the length of the MCMC chain. "    ...
                                                    + "As long as the batch size is larger than the IAC of the chain and there are significantly more than 10 "         ...
                                                    + "batches, the BatchMeans method will provide reliable estimates of the IAC." + newline + newline                  ...
                                        + "Note that in order to obtain i.i.d. samples from a multidimensional chain, " + methodName + " will use the maximum of "      ...
                                        + "IAC among all dimensions of the chain to refine the chain. Also, note that the value specified for sampleRefinementCount "   ...
                                        + "is used only when the variable sampleSize < 0, otherwise, it will be ignored. "                                              ...
                                        + "The default value is sampleRefinementMethod = '" + self.def + "'. "                                                          ...
                                        + "Note that the input values are case-insensitive and white-space characters are ignored."                                     ...
                                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, sampleRefinementMethod)

            if isempty(sampleRefinementMethod)
                self.val = self.def;
            else
                self.val = strrep(strtrim(sampleRefinementMethod), " ", "");
            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME                   = "@checkForSanity";
            sampleRefinementMethodLowerCase = lower(self.val);
            if  (index(sampleRefinementMethodLowerCase, lower(self.BATCH_MEANS)) == 0) &&             ...
                (index(sampleRefinementMethodLowerCase, lower(self.MAX_CUMSUM_AUTOCORR)) == 0)
                
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                        ...
                                + "The input requested method for the computation of the Integrated Autocorrelation Time ("     ...
                                + self.val + ") assigned to the variable sampleRefinementMethod cannot be set "                 ...
                                + "to anything other than "                                                                     ...
                                + self.BATCH_MEANS  + ". "                                                                      ...
                                ... + " or " + self.MAX_CUMSUM_AUTOCORR 
                                + "If you are not sure of the appropriate value for "                                           ...
                                + "SampleRefinementMethod, drop it from the input list. " + methodName                          ...
                                + "will automatically assign an appropriate value to it." + newline + newline                   ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end