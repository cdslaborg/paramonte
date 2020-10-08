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

function runSampler ( self          ...
                    , ndim          ...
                    , getLogFunc    ...
                    )

    self.objectName = inputname(1);

    FUNCTION_NAME = self.SUB_CLASS2_NAME + "@runParaDRAM()";

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize SpecBase variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.setupParaMonte ( ndim                      ...
                        , Constants.PMSM.MatDRAM    ...
                        , "May 23 2018"             ...
                        , self.version              ...
                        );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize ParaMCMC variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.setupParaMCMC();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize ParaDRAM variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.SpecDRAM = SpecDRAM_class( self.nd.val, self.name);     % , PD%SpecMCMC%ChainSize%def )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read variables from argument list if needed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.SpecBase.setFromInputArgs  ( ndim                                  ...
                                    , self.Err                              ...
                                    , self.spec.sampleSize                  ...
                                    , self.spec.randomSeed                  ...
                                    , self.spec.description                 ...
                                    , self.spec.outputFileName              ...
                                    , self.spec.outputDelimiter             ...
                                    , self.spec.chainFileFormat             ...
                                    , self.spec.variableNameList            ...
                                    , self.spec.restartFileFormat           ...
                                    , self.spec.outputColumnWidth           ...
                                    , self.spec.overwriteRequested          ...
                                    , self.spec.outputRealPrecision         ...
                                    , self.spec.silentModeRequested         ...
                                    , self.spec.domainLowerLimitVec         ...
                                    , self.spec.domainUpperLimitVec         ...
                                    , self.spec.parallelizationModel        ...
                                    , self.spec.progressReportPeriod        ...
                                    , self.spec.targetAcceptanceRate        ...
                                    , self.spec.maxNumDomainCheckToWarn     ...
                                    , self.spec.maxNumDomainCheckToStop     ...
                                    ) ;

    self.SpecMCMC.setFromInputArgs  ( self.SpecBase.domainLowerLimitVec.Val             ...
                                    , self.SpecBase.domainUpperLimitVec.Val             ...
                                    , self.spec.chainSize                               ...
                                    , self.spec.startPointVec                           ...
                                    , self.spec.sampleRefinementCount                   ...
                                    , self.spec.sampleRefinementMethod                  ...
                                    , self.spec.randomStartPointRequested               ...
                                    , self.spec.randomStartPointDomainLowerLimitVec     ...
                                    , self.spec.randomStartPointDomainUpperLimitVec     ...
                                    ) ;

    self.SpecDRAM.setFromInputArgs  ( self.spec.scaleFactor                             ...
                                    , self.spec.proposalModel                           ...
                                    , self.spec.proposalStartCovMat                     ...
                                    , self.spec.proposalStartCorMat                     ...
                                    , self.spec.proposalStartStdVec                     ...
                                    , self.spec.adaptiveUpdateCount                     ...
                                    , self.spec.adaptiveUpdatePeriod                    ...
                                    , self.spec.greedyAdaptationCount                   ...
                                    , self.spec.delayedRejectionCount                   ...
                                    , self.spec.burninAdaptationMeasure                 ...
                                    , self.spec.delayedRejectionScaleFactorVec          ...
                                    ) ;

    self.Err.resetEnabled   = false;
    self.Err.prefix         = self.brand;
    self.Err.newLine        = newline;
    self.Err.outputUnit     = self.LogFile.unit;

    if self.Err.occurred
        self.Err.msg        = FUNCTION_NAME + self.Err.msg;
        self.Err.abort();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now depending on the requested parallelism type, determine the process master/slave type and open output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.Image.isMaster     = false;
    if self.SpecBase.parallelizationModel.isSinglChain
        if self.Image.isFirst, self.Image.isMaster = true; end
    elseif self.SpecBase.parallelizationModel.isMultiChain
        self.Image.isMaster = true;
    else
        self.Err.msg        = FUNCTION_NAME + ": Error occurred. Unknown parallelism requested via the input variable parallelizationModel='"...
                            + self.SpecBase.parallelizationModel.val + "'.";
        self.Err.abort();
    end
    self.Image.isNotMaster  = ~self.Image.isMaster;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup log file and open output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.Chain              = ChainFileContents_class(ndim, self.SpecBase.variableNameList.Val, [], [], [], [], [], []);

    self.setupOutputFiles();


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % report variable values to the log file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if self.isFreshRun
        self.SpecBase.reportValues(self.brand, self.LogFile.unit);
        self.SpecMCMC.reportValues(self.brand, self.LogFile.unit, self.SpecBase.silentModeRequested.isFalse);
        self.SpecDRAM.reportValues(self.brand, self.LogFile.unit, self.name, self.SpecBase.silentModeRequested.isFalse);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform variable value sanity checks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.Err.msg            = "";
    self.Err.occurred       = false;

    self.SpecBase.checkForSanity (self.name, self.Err);

    self.SpecMCMC.checkForSanity (self.Err, self.name, ndim, self.SpecBase.domainLowerLimitVec.Val, self.SpecBase.domainUpperLimitVec.Val);

    self.SpecDRAM.checkForSanity (self.Err, self.name, ndim);

    if self.Err.occurred
            self.Err.msg    = FUNCTION_NAME + self.Err.msg;
            self.Err.abort();
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if self.isFreshRun
        fprintf (self.TimeFile.unit , self.TimeFile.headerFormat                ...
                                    , "NumFuncCallTotal"                        ...
                                    , "NumFuncCallAccepted"                     ...
                                    , "MeanAcceptanceRateSinceStart"            ...
                                    , "MeanAcceptanceRateSinceLastReport"       ...
                                    , "TimeElapsedSinceLastReportInSeconds"     ...
                                    , "TimeElapsedSinceStartInSeconds"          ...
                                    , "TimeRemainedToFinishInSeconds")          ...
                                    ;
        self.Chain.writeHeader ( ndim, self.ChainFile.unit, self.SpecBase.chainFileFormat.isBinary, self.ChainFile.headerFormat);
    else
        if self.Image.isMaster, fgets(self.TimeFile.unit); end   % read the header line of the time file, only by master images
        %self.Chain.getLenHeader(ndim, self.SpecBase.chainFileFormat.isBinary, self.ChainFile.headerFormat);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup the proposal distribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if self.SpecDRAM.proposalModel.isNormal || self.SpecDRAM.proposalModel.isUniform
        self.Proposal   = ParaDRAMProposalSymmetric_class(ndim, self);
    else
        self.Err.msg    = FUNCTION_NAME + ": Internal error occurred. Unsupported proposal distribution for " + self.name + ".";
        self.Err.abort();
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run ParaDRAM kernel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if self.isFreshRun && self.Image.isMaster
        self.Decor.writeDecoratedText   ( " " + newline + "Starting " + self.name + " sampling - " + DateTime_class.getNiceDateTime() + newline ...
                                        , [], [], [], [], 1, 1, self.LogFile.unit, 1 ) ;
    end

    self.runKernel(getLogFunc);

    if self.isFreshRun && self.Image.isMaster
        self.Decor.writeDecoratedText   ( " " + newline + "Exiting " + self.name + " sampling - " + DateTime_class.getNiceDateTime() + newline  ...
                                        , [], [], [], [], 1, 1, self.LogFile.unit, 1 ) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start ParaDRAM post-processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    self.Err.reset();
    self.Err.resetEnabled   = false;
    self.Err.prefix         = self.brand;
    self.Err.outputUnit     = self.LogFile.unit;
    spaces8                 = Decoration_class.TAB + Decoration_class.TAB;

    if self.Image.isMaster  % blockMasterPostProcessing

        logFileColWidthStr = num2str(max(self.SpecBase.outputRealPrecision.val, self.SpecBase.variableNameList.MaxLen.val) + 9);
        formatStr   = "%"   +       logFileColWidthStr              + "s";
        precision   = "%0." + self.SpecBase.outputRealPrecision.str + "E";
        precision2  = "%0." +                   15                  + "f";

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Total number of accepted function calls (unique samples):" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.numFuncCall.accepted" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.NumFunCall.accepted ) );
        self.Err.msg    = "This is the total number of accepted function calls (unique samples).";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Total number of accepted or rejected function calls:");
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.numFuncCall.acceptedRejected" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.NumFunCall.acceptedRejected ) );
        self.Err.msg    = "This is the total number of accepted or rejected function calls.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Total number of accepted or rejected or delayed-rejection (if any requested) function calls:");
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.numFuncCall.acceptedRejectedDelayed" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.NumFunCall.acceptedRejectedDelayed ) );
        self.Err.msg    = "This is the total number of accepted or rejected or delayed-rejection (if any requested) function calls.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "stats.chain.verbose.efficiency.meanAcceptanceRate");
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.efficiency.meanAcceptanceRate" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Chain.MeanAccRate(self.Stats.NumFunCall.accepted), precision2 ) );
        self.Err.msg    = "This is the average MCMC acceptance rate.";
        self.Err.note();

        mcmcSamplingEfficiency = self.Stats.NumFunCall.accepted / self.Stats.NumFunCall.acceptedRejected;
        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "MCMC sampling efficiency [ = acceptedFunctionCalls / acceptedPlusRejectedFunctionCalls ]:");
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.efficiency.acceptedOverAcceptedRejected" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( mcmcSamplingEfficiency, precision2 ) );
        self.Err.msg    = "This is the MCMC sampling efficiency given the accepted and rejected function calls, that is," ...
                        + "the number of accepted function calls divided by the number of (accepted + rejected) function calls.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "MCMC sampling efficiency (including delayed rejections, if any requested) [ = acceptedFunctionCalls / acceptedPlusRejectedPlusDelayedRejectionFunctionCalls ]:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.efficiency.acceptedOverAcceptedRejectedDelayed" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( (self.Stats.NumFunCall.accepted / self.Stats.NumFunCall.acceptedRejectedDelayed), precision2 ) );
        self.Err.msg    = "This is the MCMC sampling efficiency given the accepted, rejected, and delayed-rejection (if any requested) function calls, that is, " ...
                        + "the number of accepted function calls divided by the number of (accepted + rejected + delayed-rejection) function calls.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Total runtime in seconds:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.time.total" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Timer.total, precision2 ) );
        self.Err.msg    = "This is the total runtime in seconds.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Average effective time cost of each accepted function call, in seconds:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.time.perFuncCallAccepted" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Timer.total / self.Stats.NumFunCall.accepted, precision2 ) );
        self.Err.msg    = "This is the average effective time cost of each accepted function call, in seconds.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Average effective time cost of each accepted or rejected function call, in seconds:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.time.perFuncCallAcceptedRejected" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Timer.total / self.Stats.NumFunCall.acceptedRejected, precision2 ) );
        self.Err.msg    = "This is the average effective time cost of each accepted or rejected function call, in seconds.";
        self.Err.note();

        
        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Average effective time cost of each accepted or rejected or delayed-rejection function call, in seconds:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.time.perFuncCallAcceptedRejectedDelayed" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Timer.total / self.Stats.NumFunCall.acceptedRejectedDelayed ) );
        self.Err.msg    = "This is the average effective time cost of each accepted or rejected function call (including delayed-rejections, if any requested), in seconds.";
        self.Err.note();
        

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Average pure time cost of each function call, in seconds:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.time.perFuncCall" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.avgTimePerFunCalInSec ) );
        self.Err.msg    = "This is the average pure time cost of each function call, in seconds.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Number of processes (images):" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.parallelism.current.numProcess" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Image.count ) );
        self.Err.msg    = "This is the number of processes (images) used in this simulation.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Burnin location in the compact chain, based on the occurrence likelihood:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.compact.burnin.location.likelihoodBased" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.BurninLoc.compact ) );
        self.Err.msg    = "This is the burnin location in the compact chain, based on the occurrence likelihood.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Burnin location in the compact chain, based on the value of burninAdaptationMeasure:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.compact.burnin.location.adaptationBased" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.AdaptationBurninLoc.compact ) );
        self.Err.msg    = "This is the burnin location in the compact chain, based on the value of burninAdaptationMeasure simulation specification.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Burnin location in the verbose (Markov) chain, based on the occurrence likelihood:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.burnin.location.likelihoodBased" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.BurninLoc.verbose ) );
        self.Err.msg    = "This is the burnin location in the verbose (Markov) chain, based on the occurrence likelihood.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Burnin location in the verbose (Markov) chain, based on the value of burninAdaptationMeasure:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.burnin.location.adaptationBased" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.AdaptationBurninLoc.verbose ) );
        self.Err.msg    = "This is the burnin location in the verbose (Markov) chain, based on the value of burninAdaptationMeasure simulation specification.";
        self.Err.note();

        % reset BurninLoc to the maximum value

        if self.Stats.AdaptationBurninLoc.compact > self.Stats.BurninLoc.compact
            self.Stats.BurninLoc.compact = self.Stats.AdaptationBurninLoc.compact;
            self.Stats.BurninLoc.verbose = self.Stats.AdaptationBurninLoc.verbose;
        end

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Maximum-logFunc value (maximum of the user-specified objective function):" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.logFunc.max" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.LogFuncMode.val, precision2 ) );
        self.Err.msg    = "This is the maximum logFunc value (the maximum of the user-specified objective function).";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Location of the first occurrence of maximum-logFunc in the compact chain:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.compact.logFunc.max.location" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.LogFuncMode.Loc.compact ) );
        self.Err.msg    = "This is the location of the first occurrence of the maximum logFunc in the compact chain.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Location of the first occurrence of maximum-logFunc in the verbose chain:" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.logFunc.max.location" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str( self.Stats.LogFuncMode.Loc.verbose ) );
        self.Err.msg    = "This is the location of the first occurrence of the maximum logFunc in the verbose (Markov) chain.";
        self.Err.note();

        %self.Decor.write( self.LogFile.unit, 1, 0, 1, "Maximum-logFunc coordinates (mode of the user-specified objective function):" );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.logFunc.max.state" );
        self.SpecBase.variableNameList.MaxLen.val = max( self.SpecBase.variableNameList.MaxLen.val, self.SpecBase.outputColumnWidth.val );
        self.SpecBase.variableNameList.MaxLen.val = max( self.SpecBase.variableNameList.MaxLen.val, self.SpecBase.outputRealPrecision.val + 7 );
        self.SpecBase.variableNameList.MaxLen.str = num2str( self.SpecBase.variableNameList.MaxLen.val + 1 );

        fprintf( self.LogFile.unit, spaces8);
        fprintf( self.LogFile.unit, formatStr, self.SpecBase.variableNameList.Val{1:ndim} );
        fprintf( self.LogFile.unit, "\n");
        fprintf( self.LogFile.unit, spaces8);
        for i = 1 : ndim
            fprintf( self.LogFile.unit, formatStr, num2str(self.Stats.LogFuncMode.Crd(i), precision) );
        end
        self.Decor.write(self.LogFile.unit, 0, 0, [], " ");
        self.Err.msg    = "This is the coordinates, within the domain of the user-specified objective function, where the maximum NOTE: logFunc occurs.";
        self.Err.note();


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the MCMC chain's statistical properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % precision = "%0." + self.SpecBase.outputRealPrecision.val + "E";

        if self.Image.isFirst
            self.Err.msg        = "Computing the Markov chain's statistical properties...";
            self.Err.outputUnit = 1;
            self.Err.marginTop  = 3;
            self.Err.note() ;
        end

        self.Decor.writeDecoratedText   ( " " + newline + "Markov chain's statistical properties" + newline ...
                                        , [], [], [], [], 1, 1                                              ...
                                        , self.LogFile.unit                                                 ...
                                        , newline) ;

        self.Stats.Chain.count  = sum(self.Chain.Weight(self.Stats.BurninLoc.compact:self.Chain.Count.compact));

        % compute the covariance and correlation upper-triangle matrices

        self.Stats.Chain.Mean   = zeros(ndim);
        self.Stats.Chain.CovMat = zeros(ndim,ndim);
        [self.Stats.Chain.CovMat, self.Stats.Chain.Mean]    ...
                                = getWeiSamCovUppMeanTrans  ( self.Chain.Count.compact - self.Stats.BurninLoc.compact + 1                   ...
                                                            , self.Stats.Chain.count                                                        ...
                                                            , ndim                                                                          ...
                                                            , self.Chain.State(1:ndim,self.Stats.BurninLoc.compact:self.Chain.Count.compact)...
                                                            , self.Chain.Weight(self.Stats.BurninLoc.compact:self.Chain.Count.compact)      ...
                                                            ) ;

        % transpose the covariance and correlation matrices

        self.Stats.Chain.CovMat     = (self.Stats.Chain.CovMat)';

        for i = 1 : ndim
            self.Stats.Chain.CovMat(i+1:ndim,i) = self.Stats.Chain.CovMat(i,i+1:ndim);
        end

        % get the correlation matrices

        self.Stats.Chain.CorMat     = corrcov(self.Stats.Chain.CovMat);

        % compute the quantiles

        QPROB                       = QuantileProbability_class();
        self.Stats.Chain.Quantile   = zeros(QPROB.count, ndim);

        for i = 1 : ndim
            self.Stats.Chain.Quantile(1:QPROB.count,i)  = getQuantile   ( self.Chain.Count.compact - self.Stats.BurninLoc.compact + 1               ...
                                                                        , QPROB.count                                                               ...
                                                                        , QPROB.Value                                                               ...
                                                                        , self.Chain.State(i,self.Stats.BurninLoc.compact:self.Chain.Count.compact) ...
                                                                        , self.Chain.Weight(self.Stats.BurninLoc.compact:self.Chain.Count.compact)  ...
                                                                        , self.Stats.Chain.count                                                    ...
                                                                        ) ;
        end

        % report the MCMC chain statistics
        
        self.Err.marginTop          = 1;
        self.Err.outputUnit         = self.LogFile.unit;

        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.length.burninExcluded" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str(self.Stats.Chain.count) );
        self.Err.msg    = "This is the length of the verbose (Markov) Chain excluding burnin.";
        self.Err.note();

        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.avgStd" );
        fprintf(self.LogFile.unit, spaces8);
        fprintf( self.LogFile.unit, formatStr, "variableName", "Mean", "Standard Deviation", newline );
        for i = 1 : ndim
            fprintf(self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, strtrim(self.SpecBase.variableNameList.Val{i}), num2str(self.Stats.Chain.Mean(i),precision), num2str(sqrt(self.Stats.Chain.CovMat(i,i)), precision), newline );
        end
        self.Err.msg    = "This is the mean and standard deviation of the verbose (Markov) chain variables.";
        self.Err.note();

        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.covmat" );
        fprintf(self.LogFile.unit, spaces8);
        fprintf( self.LogFile.unit, formatStr, " ", self.SpecBase.variableNameList.Val{1:ndim}, newline );
        for i = 1 : ndim
            fprintf(self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, strtrim(self.SpecBase.variableNameList.Val{i}) );
            for j = 1 : ndim
                fprintf( self.LogFile.unit, formatStr, num2str(self.Stats.Chain.CovMat(i,j),precision) );
            end
            fprintf( self.LogFile.unit, formatStr, newline );
        end
        self.Err.msg    = "This is the covariance matrix of the verbose (Markov) chain.";
        self.Err.note();

        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.cormat" );
        fprintf(self.LogFile.unit, spaces8);
        fprintf( self.LogFile.unit, formatStr, " ", self.SpecBase.variableNameList.Val{1:ndim}, newline );
        for i = 1 : ndim
            fprintf(self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, strtrim(self.SpecBase.variableNameList.Val{i}) );
            for j = 1 : ndim
                fprintf( self.LogFile.unit, formatStr, num2str(self.Stats.Chain.CorMat(i,j),precision) );
            end
            fprintf( self.LogFile.unit, formatStr, newline );
        end
        self.Err.msg    = "This is the correlation matrix of the verbose (Markov) chain.";
        self.Err.note();

        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.verbose.quantile" );
        fprintf(self.LogFile.unit, spaces8);
        fprintf( self.LogFile.unit, formatStr, "Quantile", self.SpecBase.variableNameList.Val{1:ndim}, newline );
        for iq = 1 : QPROB.count
            fprintf(self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, num2str(QPROB.Name(iq),precision) );
            for i = 1 : ndim
            fprintf( self.LogFile.unit, formatStr, num2str(self.Stats.Chain.Quantile(iq,i),precision) );
            end
            fprintf( self.LogFile.unit, formatStr, newline );
        end
        self.Err.msg    = "This are the quantiles table of the variables of the verbose (Markov) chain.";
        self.Err.note();

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate the i.i.d. sample statistics and output file (if requested)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % report refined sample statistics, and generate output refined sample if requested.

        % if self.Image.isFirst
        %     self.Err.msg    = "Computing the final decorrelated sample size...";
        %     self.Err.note();
        % end

       [self.Err, self.Chain]   = self.RefinedChain.get ( self.Chain                                ...
                                                        , self.Stats.BurninLoc.compact              ...
                                                        , []                                        ...
                                                        , self.SpecMCMC.sampleRefinementCount.val   ...
                                                        , self.SpecMCMC.sampleRefinementMethod.val  ...
                                                        ) ;
        self.Err.resetEnabled   = false;
        self.Err.prefix         = self.brand;
        self.Err.outputUnit     = self.LogFile.unit;
        
        if self.Err.occurred
            self.Err.msg        = FUNCTION_NAME + self.Err.msg;
            self.Err.abort();
        end

        % compute the maximum integrated autocorrelation times for each variable

        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.iac" );
        fprintf(self.LogFile.unit, spaces8);
        fprintf(self.LogFile.unit, formatStr, "RefinementStage", "SampleSize", "IAC(SampleLogFunc)" );
        for i = 1 : ndim
            fprintf( self.LogFile.unit, formatStr, "IAC(" + strtrim(self.SpecBase.variableNameList.Val{i}) + ")" );
        end
        self.Decor.write( self.LogFile.unit, 0, 0, [], " " );
        for i = 0 : self.RefinedChain.numRefinement
            fprintf(self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, num2str(i), num2str(self.RefinedChain.Count(i+1).verbose) );
            for j = 1 : ndim+1
                fprintf( self.LogFile.unit, formatStr, num2str(self.RefinedChain.IAC(j,i+1),precision) );
            end
            fprintf( self.LogFile.unit, formatStr, newline );
        end
        self.Err.msg    = "This is the table of the Integrated Autocorrelation (IAC) of individual variables in the verbose (Markov) chain, at increasing stages of chain refinements.";
        self.Err.note();

        % report the final Effective Sample Size (ESS) based on IAC

        effectiveSampleSize = sum( self.RefinedChain.Weight(1:self.RefinedChain.Count(self.RefinedChain.numRefinement+1).compact) );
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.ess");
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str(effectiveSampleSize) );
        self.Err.msg    = "This is the estimated Effective (decorrelated) Sample Size (ESS) of the final refined chain.";
        self.Err.note();
        
        self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.efficiency.essOverAccepted" );
        self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str(effectiveSampleSize / self.Stats.NumFunCall.accepted,precision2) );
        self.Err.msg    = "This is the effective MCMC sampling efficiency given the accepted function calls, that is, the final " ...
                        + "refined effective sample size (ESS) divided by the number of accepted function calls.";
        self.Err.note();

        % if self.SpecDRAM.delayedRejectionCount.val == 0
            self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.efficiency.essOverAcceptedRejected" );
            self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str(effectiveSampleSize / self.Stats.NumFunCall.acceptedRejected,precision2) );
            self.Err.msg    = "This is the effective MCMC sampling efficiency given the accepted and rejected function calls, that " ...
                            + " is, the final refined effective sample size (ESS) divided by the number of (accepted + rejected) " ...
                            + "function calls.";
            self.Err.note();
        % else
            self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.efficiency.essOverAcceptedRejectedDelayed" );
            self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str(effectiveSampleSize / self.Stats.NumFunCall.acceptedRejectedDelayed,precision2) );
            self.Err.msg    = "This is the effective MCMC sampling efficiency given the accepted, rejected, and delayed-rejection " ...
                            + "(if any requested) function calls, that is, the final refined effective sample size (ESS) divided by " ...
                            + "the number of (accepted + rejected + delayed-rejection) function calls.";
            self.Err.note();
        % end

        % generate output refined sample if requested

        if self.SpecBase.sampleSize.val == 0    % blockSampleFileGeneration
            self.Err.msg    = "Skipping decorrelated sample generation and output file, as requested by user...";
            self.Err.note();

        else % blockSampleFileGeneration

            %*******************************************************************************************************************
            % sample file generation report
            %*******************************************************************************************************************

            % report to the report-file(s)
            
            self.Err.msg    = "Generating the output " + self.SampleFile.suffix + " file:" + newline + strrep(self.SampleFile.Path.original, '\', '\\');
            self.Err.note();

            if self.Image.isFirst

                % print the message for the generating the output sample file on the first image
                
                self.Err.outputUnit = 1;
                
                self.Err.msg        = "Generating the output " + self.SampleFile.suffix + " file:";
                self.Err.marginBot  = 0;
                self.Err.note();
                
                self.Err.msg        = strrep(self.SampleFile.Path.original, '\', '\\');
                self.Err.marginTop  = 0;
                self.Err.marginBot  = 0;
                self.Err.note();

            end

            %*******************************************************************************************************************
            % begin sample file generation
            %*******************************************************************************************************************

            if self.SpecBase.sampleSize.val ~= -1

                if self.SpecBase.sampleSize.val < 0, self.SpecBase.sampleSize.val = abs(self.SpecBase.sampleSize.val) * self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose; end
                
                self.Err.outputUnit = self.LogFile.unit;
                if self.SpecBase.sampleSize.val < self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose
                    self.Err.msg    = "The user-requested sample size (" + num2str(self.SpecBase.sampleSize.val) + ") "         ...
                                    + "is smaller than the optimal i.i.d. sample size "                                         ...
                                    + "(" + num2str(self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose) + "). " ...
                                    + "The output sample contains i.i.d. samples, however, the sample-size "                    ...
                                    + "could have been larger if it had been set to the optimal size. "                         ...
                                    + "To get the optimal size in the future runs, set sampleSize = -1, or drop "               ...
                                    + "it from the input list."                                                                 ...
                                    ;
                    self.Err.warn() ;
                elseif self.SpecBase.sampleSize.val > self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose
                    self.Err.msg    = "The user-requested sample size (" + num2str(self.SpecBase.sampleSize.val) + ") "         ...
                                    + "is larger than the optimal i.i.d. sample size "                                          ...
                                    + "(" + num2str(self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose) + "). " ...
                                    + "The resulting sample likely contains duplicates and is not independently "               ...
                                    + "and identically distributed (i.i.d.)." + newline + "To get the optimal "                 ...
                                    + "size in the future runs, set sampleSize = -1, or drop "                                  ...
                                    + "it from the input list."                                                                 ...
                                    ;
                    self.Err.warn() ;
                else
                    self.Err.msg    = "How lucky that could be! "                                                               ...
                                    + "The user-requested sample size (" + num2str(self.SpecBase.sampleSize.val) + ") "         ...
                                    + "is equal to the optimal i.i.d. sample size determined by " + self.name + "."             ...
                                    ;
                    self.Err.warn() ;
                end

                % regenerate the refined sample, this time with the user-requested sample size.

                [self.Err, self.Chain] = self.RefinedChain.get  ( self.Chain                                ...
                                                                , self.Stats.BurninLoc.compact              ...
                                                                , self.SpecBase.sampleSize.val              ...
                                                                , self.SpecMCMC.sampleRefinementCount.val   ...
                                                                , self.SpecMCMC.sampleRefinementMethod.val  ...
                                                                ) ;
                if self.Err.occurred
                    self.Err.msg    = FUNCTION_NAME + self.Err.msg;
                    self.Err.abort();
                end

            end

            % open the output sample file

            [self.SampleFile.unit, self.Err.msg] = fopen(self.SampleFile.Path.original, "W");
            if self.Err.msg
                self.Err.msg    = FUNCTION_NAME + ": Error occurred while opening " + self.name + " " + self.SampleFile.suffix + " file='" + self.SampleFile.Path.original + "'.";
                self.Err.abort();
            end

            % write to the output sample file

            self.RefinedChain.write(self.SampleFile.unit, self.SampleFile.headerFormat, self.SampleFile.format);
            fclose(self.SampleFile.unit);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute the refined sample's statistical properties
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if self.Image.isFirst
                self.Err.msg        = "Computing the output sample's statistical properties...";
                self.Err.marginTop  = 2;
                self.Err.outputUnit = 1;
                self.Err.note();
            end

            self.Decor.writeDecoratedText   ( " " + newline + "The statistical properties of the final refined sample" + newline    ...
                                            , [], [], [], [], 1, 1                                                  ...
                                            , self.LogFile.unit                                                     ...
                                            , newline) ;

            self.Stats.Sample.count     = self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose;   % xxx

            % compute the covariance and correlation upper-triangle matrices

            self.Stats.Sample.Mean      = zeros(1, ndim);
            self.Stats.Sample.CovMat    = zeros(ndim,ndim);
            [self.Stats.Sample.CovMat, self.Stats.Sample.Mean]  = getWeiSamCovUppMeanTrans  ( self.RefinedChain.Count(self.RefinedChain.numRefinement+1).compact                ...
                                                                , self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose                                            ...
                                                                , ndim                                                                                                          ...
                                                                , self.RefinedChain.LogFuncState(2:ndim+1,1:self.RefinedChain.Count(self.RefinedChain.numRefinement+1).compact) ...
                                                                , self.RefinedChain.Weight(1:self.RefinedChain.Count(self.RefinedChain.numRefinement+1).compact)                ...
                                                                ) ;

            % transpose the covariance matrix

            self.Stats.Sample.CovMat    = (self.Stats.Sample.CovMat)';

            for i = 1 : ndim
                self.Stats.Sample.CovMat(i+1:ndim,i) = self.Stats.Sample.CovMat(i,i+1:ndim);
            end

            % get the correlation matrices

            self.Stats.Sample.CorMat    = corrcov(self.Stats.Sample.CovMat);

            % compute the quantiles

            self.Stats.Sample.Quantile = zeros(QPROB.count, ndim);

            for i = 1 : ndim
                self.Stats.Sample.Quantile(1:QPROB.count,i) = getQuantile   ( self.RefinedChain.Count(self.RefinedChain.numRefinement+1).compact                                        ...
                                                                            , QPROB.count                                                                                               ...
                                                                            , QPROB.Value                                                                                               ...
                                                                            , self.RefinedChain.LogFuncState(i+1,1:self.RefinedChain.Count(self.RefinedChain.numRefinement+1).compact)  ...
                                                                            , self.RefinedChain.Weight(1:self.RefinedChain.Count(self.RefinedChain.numRefinement+1).compact)            ...
                                                                            , self.RefinedChain.Count(self.RefinedChain.numRefinement+1).verbose                                        ...
                                                                            ) ;
            end

            % report the MCMC chain statistics
            
            self.Err.outputUnit     = self.LogFile.unit;
            self.Err.marginTop      = 1;
            self.Err.marginBot      = 1;

            self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.length" );
            self.Decor.write( self.LogFile.unit, 0, 0, 1, spaces8 + num2str(self.Stats.Sample.count) );
            self.Err.msg    = "This is the final output refined sample size.";
            self.Err.note();

            self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.avgStd" );
            fprintf( self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, "variablename", "Mean", "Standard Deviation", newline );
            for i = 1 : ndim
                fprintf( self.LogFile.unit, spaces8);
                fprintf( self.LogFile.unit, formatStr, strtrim(self.SpecBase.variableNameList.Val{i}), num2str(self.Stats.Sample.Mean(i),precision), num2str(sqrt(self.Stats.Sample.CovMat(i,i)), precision), newline );
            end
            self.Err.msg    = "This is the Mean and standard deviation table of the final output refined sample.";
            self.Err.note();

            self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.covmat" );
            fprintf( self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, "", self.SpecBase.variableNameList.Val{1:ndim}, newline );
            for i = 1 : ndim
                fprintf( self.LogFile.unit, spaces8);
                fprintf( self.LogFile.unit, formatStr, strtrim(self.SpecBase.variableNameList.Val{i}) );
                for j = 1 : ndim
                    fprintf( self.LogFile.unit, formatStr, num2str(self.Stats.Sample.CovMat(i,j),precision) );
                end
                fprintf( self.LogFile.unit, formatStr, newline );
            end
            self.Err.msg    = "This is the covariance matrix of the final output refined sample.";
            self.Err.note();

            self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.cormat" );
            fprintf( self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, "", self.SpecBase.variableNameList.Val{1:ndim}, newline );
            for i = 1 : ndim
                fprintf( self.LogFile.unit, spaces8);
                fprintf( self.LogFile.unit, formatStr, strtrim(self.SpecBase.variableNameList.Val{i}) );
                for j = 1 : ndim
                    fprintf( self.LogFile.unit, formatStr, num2str(self.Stats.Sample.CorMat(i,j),precision) );
                end
                fprintf( self.LogFile.unit, formatStr, newline );
            end
            self.Err.msg    = "This is the correlation matrix of the final output refined sample.";
            self.Err.note();

            self.Decor.write( self.LogFile.unit, 1, 1, 1, "stats.chain.refined.quantile" );
            fprintf( self.LogFile.unit, spaces8);
            fprintf( self.LogFile.unit, formatStr, "Quantile", self.SpecBase.variableNameList.Val{1:ndim}, newline );
            for iq = 1 : QPROB.count
                fprintf( self.LogFile.unit, spaces8);
                fprintf( self.LogFile.unit, formatStr, num2str(QPROB.Name(iq),precision) );
                for i = 1 : ndim
                fprintf( self.LogFile.unit, formatStr, num2str(self.Stats.Sample.Quantile(iq,i),precision) );
                end
                fprintf( self.LogFile.unit, formatStr, newline );
            end
            self.Err.msg    = "This is the quantiles table of the variables of the final output refined sample.";
            self.Err.note();


        end % blockSampleFileGeneration

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % End of generating the i.i.d. sample statistics and output file (if requested)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        self.Decor.write( self.LogFile.unit, 2, 0, 0, "" );

        % Mission accomplished.
        
        self.Err.msg            = "Mission Accomplished.";

        self.Err.outputUnit     = self.LogFile.unit;
        self.Err.note();

        if 1 ~= self.LogFile.unit && self.Image.isFirst
            self.Err.msg        = "Mission Accomplished.";
            self.Err.outputUnit = 1;
            self.Err.note();
            self.Decor.write( 1, 1, 1, [], []);
        end

        fclose(self.TimeFile.unit);
        fclose(self.LogFile.unit);

    end % blockMasterPostProcessing

    self.Err.msg    = "To read the generated output files, try the following:" + newline ...
                    + newline ...
                    + "    " + self.objectName + ".readReport()      %% to read the summary report from the output report file." + newline ...
                    + "    " + self.objectName + ".readSample()      %% to read the final i.i.d. sample from the output sample file." + newline ...
                    + "    " + self.objectName + ".readChain()       %% to read the uniquely-accepted points from the output chain file." + newline ...
                    + "    " + self.objectName + ".readMarkovChain() %% to read the Markov Chain. NOT recommended for extremely-large chains." + newline ...
                    + "    " + self.objectName + ".readRestart()     %% to read the contents of an ASCII-format output restart file." + newline ...
                    + newline ...
                    + "For more information and examples on the usage, visit:" + newline ...
                    + newline ...
                    + "    " + href(self.website.home.url);
    self.Err.note();

end % function runSampler
