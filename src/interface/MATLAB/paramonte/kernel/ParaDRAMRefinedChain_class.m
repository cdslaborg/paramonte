%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ParaDRAMRefinedChain_class < handle

    properties (Constant)
        CLASS_NAME = "@ParaDRAMRefinedChain_mod"
    end

    properties
        ndim            = 0             % number of sampling variables
        numRefinement   = 0             % number of refinements, zero if sample size is prescribed by the user
        Count           = Count_class() % compact and verbose counts
        IAC             = []            % size of (ndim,0:numRefinement): The Integrated AutoCorrelation Time at each refinement stage
        LogFuncState    = []            % size of (ndim,Count%compact): LogFuncState is LogFunc + Variables
        Weight          = []            % size of (Count%compact): Weight of each state
        ColHeader       = [""]          % refined sample column headers
        Err             = Err_class()
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function [Err, CFC] =   get ( self                      ...
                                    , CFC                       ...
                                    , burninLoc                 ...
                                    , refinedChainSize          ...
                                    , sampleRefinementCount     ...
                                    , sampleRefinementMethod    ...
                                    )
            FUNCTION_NAME = self.CLASS_NAME + "@getRefinedChain()";

            Err             = Err_class();
            Err.occurred    = false;
            Method          = Method_class();

            % if the size of the refined sample is given as input, then generate the requested sample straight

            refinedChainSizeIsPresent = ~isempty(refinedChainSize);
            if refinedChainSizeIsPresent    % ignore sampleRefinementCount, even if it is given by the user
                maxRefinementCount  = 1;
            else
                maxRefinementCount  = 20;   % this is a temporary maximum value, to be increased later if needed
                if ~isempty(sampleRefinementCount), maxRefinementCount = sampleRefinementCount; end
            end

            % this is to avoid memory overflow due to extremely large maxRefinementCount requested by the user

            maxRefinementCurrent    = min(2, maxRefinementCount);

            % compute ndim and the initial chain size

            self.numRefinement      = 0;

            if CFC.ndim == 0
                self.ndim   = length(CFC.State(:,1));
                CFC.ndim    = self.ndim;
            else
                self.ndim   = CFC.ndim;
            end

            % allocate components

            self.IAC                                = zeros(self.ndim+1, maxRefinementCurrent+1);
            self.Count(1:maxRefinementCurrent+1)    = Count_class();

            if CFC.Count.compact == 0
                self.Count(self.numRefinement+1).compact    = length(CFC.State(1,:));
                CFC.Count.compact                           = self.Count(self.numRefinement+1).compact;
            else
                self.Count(self.numRefinement+1).compact    = CFC.Count.compact;
            end
            if CFC.Count.verbose == 0, CFC.Count.verbose    = sum(CFC.Weight(1:CFC.Count.compact)); end
            self.Count(self.numRefinement+1).verbose        = CFC.Count.verbose;

            BATCH_MEANS         = SpecMCMC_SampleRefinementMethod_class.BATCH_MEANS;
            MAX_CUMSUM_AUTOCORR = SpecMCMC_SampleRefinementMethod_class.MAX_CUMSUM_AUTOCORR;
            if ~isempty(sampleRefinementMethod)
                sampleRefinementMethodLowerCase = lower(sampleRefinementMethod);
                if      index(sampleRefinementMethodLowerCase, lower(BATCH_MEANS)) > 0
                    Method.isBatchMeans         = true;
                elseif  index(sampleRefinementMethodLowerCase, lower(MAX_CUMSUM_AUTOCORR)) > 0
                    Method.isMaxCumSumAutoCorr  = true;
                else
                    Err.occurred    = true;
                    Err.msg         = FUNCTION_NAME + ": Unknown unsupported IAC computation method name: " + sampleRefinementMethod;
                    return
                end
                Method.isViaCompactChain = index(sampleRefinementMethodLowerCase, "compact") > 0;
                Method.isViaVerboseChain = index(sampleRefinementMethodLowerCase, "verbose") > 0;
                if ~(Method.isViaCompactChain || Method.isViaVerboseChain)
                    Method.isViaCompactChain = true;
                    Method.isViaVerboseChain = true;
                end
            else
                Method.isBatchMeans = true; % this is the default method
            end

            % assign the column headers

            NUM_DEF_COL = ChainFileContents_class.NUM_DEF_COL;
            if ~isempty(CFC.ColHeader)
                for i = 0 : self.ndim
                    self.ColHeader(i+1) = CFC.ColHeader(i+NUM_DEF_COL);
                end
            else
                Err.occurred    = true;
                Err.msg         = FUNCTION_NAME + ": Internal error occurred. CFC.ColHeader is empty.";
                return
            end

            if ~isempty(burninLoc)
                burninLocDefault = burninLoc;
            else
                burninLocDefault = CFC.BurninLoc(CFC.Count.compact);
            end

            self.Count(self.numRefinement+1).compact = CFC.Count.compact - burninLocDefault + 1;
            self.LogFuncState                                                               = zeros(self.ndim+1, self.Count(self.numRefinement+1).compact);
            self.LogFuncState(1            , 1:self.Count(self.numRefinement+1).compact)    = CFC.LogFunc(burninLocDefault:CFC.Count.compact);
            self.LogFuncState(2:self.ndim+1, 1:self.Count(self.numRefinement+1).compact)    = CFC.State  (1:self.ndim, burninLocDefault:CFC.Count.compact);

            % check if there are more than 1 sample points in the burnin-subtracted CFC

            if self.Count(self.numRefinement+1).compact == 0
                Err.occurred    = true;
                Err.msg         = FUNCTION_NAME + ": The size of the refined sample is zero.";
                return
            elseif self.Count(self.numRefinement+1).compact == 1
                self.Weight = zeros(1, self.Count(self.numRefinement+1).compact);
                if refinedChainSizeIsPresent
                    self.Weight(:) = refinedChainSize;
                else
                    self.Weight(:) = 1;
                end
                self.IAC    = zeros(1:self.ndim+1,1);
                self.Count  = Count_class();
                self.Count(self.numRefinement+1).verbose = sum(self.Weight);
                return
            end

            % self.Weight = zeros( 1,length(CFC.Weight(burninLocDefault:CFC.Count.compact)) );    % xxx
            % self.Weight = CFC.Weight;
            self.Weight = CFC.Weight(burninLocDefault:CFC.Count.compact);   % Weight is intentionally separately assigned from State here

            while true  % loopRefinement

                countCompact = self.Count(self.numRefinement+1).compact;

                % set the sample weight

                if Method.isViaCompactChain
                    SampleWeight    = [];
                elseif Method.isViaVerboseChain
                    SampleWeight    = self.Weight;
                end

                % obtain the IAC for each individual variable

                for i = 0 : self.ndim
                    DumVec = self.LogFuncState(i+1,1:self.Count(self.numRefinement+1).compact);
                    if Method.isBatchMeans
                        self.IAC(i+1,self.numRefinement+1) = getBatchMeansIAC   ( countCompact      ...
                                                                                , DumVec            ...
                                                                                , SampleWeight      ...
                                                                                , []                ...
                                                                                ) ;
                    elseif Method.isMaxCumSumAutoCorr
                        self.IAC(i+1,self.numRefinement+1) = getMaxCumSumIAC    ( countCompact      ...
                                                                                , DumVec            ...
                                                                                , SampleWeight      ...
                                                                                ) ;
                    end
                end

                if refinedChainSizeIsPresent
                    integratedAutoCorrTime = sum(self.Weight) / refinedChainSize;
                else
                    integratedAutoCorrTime = max( self.IAC(1:self.ndim+1,self.numRefinement+1) );
                end

                % so far, we have computed the max IAC of the sample, but no refinement. Refine the sample only if needed.

                maxRefinementCountIsReached = self.numRefinement == maxRefinementCount;

                if integratedAutoCorrTime < 2 || maxRefinementCountIsReached
                    if Method.isViaCompactChain && Method.isViaVerboseChain
                        if maxRefinementCountIsReached, maxRefinementCount = maxRefinementCount * 2; end
                        Method.isViaCompactChain = false;
                        continue
                    end
                    break
                end

                % generate the refined sample, dump it in CFC, then put it back into RefinedChain to start over again

                self.numRefinement = self.numRefinement + 1;

                % reallocate to bigger array if nedded

                if self.numRefinement > maxRefinementCurrent
                    DumIAC                                          = self.IAC;
                    DumCount                                        = self.Count;
                    maxRefinementCurrent                            = min( maxRefinementCurrent*2, maxRefinementCount );
                    self.IAC                                        = zeros( self.ndim+1, maxRefinementCurrent+1);
                    self.Count(1:maxRefinementCurrent+1)            = Count_class();
                    self.IAC(1:self.ndim+1, 1:self.numRefinement)   = DumIAC;
                    self.Count(1:self.numRefinement)                = DumCount;
                end

                if integratedAutoCorrTime < 2, continue; end    % no need for refinement. should happen only when transitioning from compact to verbose

                [DumCFC.State, DumCFC.Weight, self.Count(self.numRefinement+1)] = ...
                                        ParaDRAMRefinedChain_class.refineWeightedSample ( self.ndim                             ...
                                                                                        , countCompact                          ...
                                                                                        , integratedAutoCorrTime                ...
                                                                                        , self.LogFuncState                     ...
                                                                                        , self.Weight                           ...
                                                                                        , refinedChainSize                      ...
                                                                                        ) ;
                self.Weight        = DumCFC.Weight;
                self.LogFuncState  = DumCFC.State;

            end % loopRefinement

        end % get

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function write  ( self                      ...
                        , sampleFileUnit            ...
                        , sampleFileHeaderFormat    ...
                        , sampleFileContentsFormat  ...
                        )

            fprintf( sampleFileUnit, sampleFileHeaderFormat,    self.ColHeader(:) );

            for isample = 1 : self.Count(self.numRefinement+1).compact
                for iweight = 1 : self.Weight(isample)
                    fprintf( sampleFileUnit, sampleFileContentsFormat,  self.LogFuncState(:,isample) );
                end
            end
            
            % for isample = 1 : self.Count(self.numRefinement+1).verbose
                % fprintf( sampleFileUnit, sampleFileContentsFormat,  self.LogFuncState(:,isample) );
            % end

            % for isample = 1, RefinedChain.Count(RefinedChain.numRefinement).compact
                % for iweight = 1, RefinedChain.Weight(isample)
                    % write(sampleFileUnit,sampleFileContentsFormat) RefinedChain.LogFuncState(0:RefinedChain.ndim,isample)
                % end
            % end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Static)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function [RefinedChain, RefinedWeight, PointCount] = refineWeightedSample(nd, np, skip, Sample, Weight, refinedChainSize)
            UpdatedWeight       = ParaDRAMRefinedChain_class.getRefinedWeight(np,Weight,skip,refinedChainSize);
            npRefined           = sum(UpdatedWeight > 0);
            RefinedChain        = zeros(nd+1,npRefined);
            RefinedWeight       = zeros(1,npRefined);
            PointCount          = Count_class();
            ipRefined           = 0;
            PointCount.verbose  = 0;
            for ip = 1 : np
                if UpdatedWeight(ip) > 0
                    ipRefined = ipRefined + 1;
                    RefinedChain(1:nd+1,ipRefined) = Sample(1:nd+1,ip);
                    RefinedWeight(ipRefined) = UpdatedWeight(ip);
                    PointCount.verbose = PointCount.verbose + RefinedWeight(ipRefined);
                end
            end
            PointCount.compact  = npRefined;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function RefinedWeight = getRefinedWeight(np,Weight,skip,refinedChainSize)
            RefinedWeight               = zeros(1, np);
            refinedChainSizeIsPresent   = ~isempty(refinedChainSize);
            if refinedChainSizeIsPresent, refinedChainSizeCounter = 0; end
            CumSumWeight    = cumsum(Weight);
            sumSkips            = skip;
            offset              = 1;
            ip                  = offset;
            RefinedWeight(:)    = 0;
            while true % loopOverAllSample
                while true % loopOverCurrentSample

                    if sumSkips > CumSumWeight(ip)
                        break % loopOverCurrentSample
                    elseif refinedChainSizeIsPresent
                        if refinedChainSizeCounter == refinedChainSize, break; end  % loopOverAllSample
                        refinedChainSizeCounter = refinedChainSizeCounter + 1;
                    end
                    RefinedWeight(ip)   = RefinedWeight(ip) + 1;
                    sumSkips            = sumSkips + skip;
                    continue % loopOverCurrentSample
                end %loopOverCurrentSample
                if ip == np
                    if refinedChainSizeIsPresent
                        if refinedChainSizeCounter < refinedChainSize
                            offset = offset + 1;
                            if offset == np, offset = 1; end
                            ip          = offset;
                            sumSkips    = skip;
                            if offset ~= 1, sumSkips = sumSkips + CumSumWeight(ip-1); end
                            continue % loopOverAllSample
                        end
                    end
                    break % loopOverAllSample
                end
                ip = ip + 1;
            end % loopOverAllSample
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end