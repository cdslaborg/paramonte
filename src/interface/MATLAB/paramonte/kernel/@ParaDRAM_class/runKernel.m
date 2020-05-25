function runKernel  ( self          ...
                    , getLogFunc    ...
                    )

    FUNCTION_NAME           = self.SUB_CLASS2_NAME + "@runKernel()";
    CHAIN_RESTART_OFFSET    = 2;

    global timeElapsedUntilLastReportInSeconds SumAccRateSinceStart sumAccRateLastReport inverseProgressReportPeriod numFunCallAcceptedRejectedLastReport

    global co_proposalFound_samplerUpdateOccurred   % merging these scalars would reduce the MPI communication
                                                    % overhead cost: co_proposalFound, co_samplerUpdateOccurred,
                                                    % co_counterDRS, 0 means false, 1 means true

    SumAccRateSinceStart                        = SumAccRateSinceStart_class();

    acceptedRejectedDelayedUnusedRestartMode    = 0;                                            % used to compute more accurate timings in the restart mode
    self.Stats.avgTimePerFunCalInSec            = 0.0;
    numFunCallAcceptedRejectedLastReport        = 0;
    timeElapsedUntilLastReportInSeconds         = 0.0;
    inverseProgressReportPeriod                 = 1.0 / self.SpecBase.progressReportPeriod.val; % this remains a constant except for the last the last report of the simulation
    sumAccRateLastReport                        = 0.0;
    nd                                          = self.nd.val;

    co_LogFuncState = zeros(nd+1, self.SpecDRAM.delayedRejectionCount.val+2);
    co_AccRate      = zeros(1   , self.SpecDRAM.delayedRejectionCount.val+2);

    co_AccRate(1)                                           = 0;  % the real-value counterDRS, indicating the initial delayed rejection stage at which the first point is sampled
    co_AccRate(2)                                           = 1;  % initial acceptance rate for the first zeroth DR stage.
    co_AccRate(3:self.SpecDRAM.delayedRejectionCount.val+2) = 0;  % indicates the very first proposal acceptance on image 1

    delayedRejectionRequested               = self.SpecDRAM.delayedRejectionCount.val > 0;
    noDelayedRejectionRequested             = ~delayedRejectionRequested;

    if delayedRejectionRequested
        self.Stats.NumFunCall.acceptedRejectedDelayed   = 0;    % Markov Chain counter
        SumAccRateSinceStart.acceptedRejectedDelayed    = 0.0;  % sum of acceptance rate
    end

    SumAccRateSinceStart.acceptedRejected   = 0.0;  % sum of acceptance rate
    self.Stats.NumFunCall.acceptedRejected  = 0;    % Markov Chain counter
    counterAUC                              = 0;    % counter for padaptiveUpdateCount.
    counterPRP                              = 0;    % counter for progressReportPeriod.
    counterAUP                              = 0;    % counter for adaptiveUpdatePeriod.
    self.Stats.NumFunCall.accepted          = 0;    % Markov Chain acceptance counter.
    samplerUpdateSucceeded                  = true; % needed to set up lastStateWeight and numFunCallAcceptedLastAdaptation for the first accepted proposal
    chainAdaptationMeasure                  = 0.0;  % needed for the first output
    numFunCallAcceptedLastAdaptation        = 0;
    lastStateWeight                         = -intmax;

    self.Timer.tic;
    
    self.Err.prefix     = self.brand;
    self.Err.outputUnit = self.LogFile.unit;

    if self.isFreshRun  % blockDryRunSetup

        self.Chain.ProcessID    = zeros (1  , self.SpecMCMC.chainSize.val);
        self.Chain.DelRejStage  = zeros (1  , self.SpecMCMC.chainSize.val);
        self.Chain.MeanAccRate  = zeros (1  , self.SpecMCMC.chainSize.val);
        self.Chain.Adaptation   = zeros (1  , self.SpecMCMC.chainSize.val);
        self.Chain.BurninLoc    = zeros (1  , self.SpecMCMC.chainSize.val);
        self.Chain.Weight       = zeros (1  , self.SpecMCMC.chainSize.val);
        self.Chain.State        = zeros (nd , self.SpecMCMC.chainSize.val);
        self.Chain.LogFunc      = zeros (1  , self.SpecMCMC.chainSize.val);

    else % blockDryRunSetup

        % load the existing Chain file into PD%Chain components

        self.Err = self.Chain.get   ( self.ChainFile.Path.original      ...
                                    , self.SpecBase.chainFileFormat.val ...
                                    , self.SpecMCMC.chainSize.val       ...
                                    , self.Chain.lenHeader              ...
                                    , nd                                ...
                                    , self.SpecBase.outputDelimiter.val ...
                                    , []                                ...
                                    ) ;

        if self.Err.occurred
            self.Err.msg    = FUNCTION_NAME + self.Err.msg;
            self.Err.abort();
        end

        if self.Chain.Count.compact <= CHAIN_RESTART_OFFSET
            self.isFreshRun = true;
            self.isDryRun   = ~self.isFreshRun;
        end

        if self.Image.isMaster

            % set up the chain file

            % create a copy of the chain file, just for the sake of not losing the simulation results

            RFN         = RandomFileName_class("", self.ChainFile.Path.original + ".temporary_restart_copy", "");
            self.Err    = copyFile(self.ChainFile.Path.original, RFN.path, self.OS.isWindows);

            if self.Err.occurred
                self.Err.msg    = FUNCTION_NAME + self.Err.msg;
                self.Err.abort();
            end

            % reopen the chain file to resume the simulation

            [self.ChainFile.unit, self.Err.msg] = fopen(self.ChainFile.Path.original, "w");
            if self.Err.msg
                self.Err.msg    = FUNCTION_NAME + ": Error occurred while opening " + self.name + self.ChainFile.suffix + " file='" + self.ChainFile.Path.original + "'.";
                self.Err.abort();
            end

            % rewrite the chain file

            self.Chain.writeChainFile   ( nd                                                ...
                                        , 1                                                 ...
                                        , self.Chain.Count.compact - CHAIN_RESTART_OFFSET   ...
                                        , self.ChainFile.unit                               ...
                                        , self.SpecBase.chainFileFormat.val                 ...
                                        , self.ChainFile.format                             ...
                                        , self.ChainFile.headerFormat                       ...
                                        ) ;

%************************************************************************************************************
% Restart mode old chain ends here
fprintf(self.ChainFile.unit, "-----------------Old chain with less two points till here.---------------\n");
%************************************************************************************************************

            % remove the temporary copy of the chain file

            self.Err = removeFile(RFN.path, self.OS.isWindows);
            if self.Err.occurred
                self.Err.msg    = FUNCTION_NAME + self.Err.msg;
                self.Err.abort();
            end

        end

    end % blockDryRunSetup

    if self.isFreshRun  % this must be done separately from the above blockDryRunSetup
        self.Chain.BurninLoc(1)  = 1;
        co_LogFuncState(2:nd+1,2)   = self.SpecMCMC.startPointVec.Val;              % proposal state
        self.Timer.toc();
        co_LogFuncState(1,2)        = getLogFunc(co_LogFuncState(2:nd+1,2));    % proposal logFunc
        self.Timer.toc();
        self.Stats.avgTimePerFunCalInSec = self.Stats.avgTimePerFunCalInSec + self.Timer.delta;
    else
        co_LogFuncState(2:nd+1,2)   = self.Chain.State(1:nd,1); % proposal logFunc
        co_LogFuncState(1,2)        = self.Chain.LogFunc(1);    % proposal logFunc
    end

    co_LogFuncState(1:nd+1,1)           = co_LogFuncState(1:nd+1,2);  % set current logFunc and State equal to the first proposal
    self.Stats.LogFuncMode.val          = -intmax;
    self.Stats.LogFuncMode.Loc.compact  = 0;

    if self.Image.isFirst
        fprintf( 1, "            Accepted/Total Func. Call   Dynamic/Overall Acc. Rate   Elapsed/Remained Time [s]\n"   );
        fprintf( 1, "            =========================   =========================   =========================\n"   );
        fprintf( 1, "                                                                                             "     );
    end

    imageID = 1;    % needed even in the case of serial run to assign a proper value to self.Chain.ProcessID

    %***********************************************************************************************************************
    %***********************************************************************************************************************
    %******************************************** start of loopMarkovChain *************************************************

    while true  % loopMarkovChain

        co_proposalFound_samplerUpdateOccurred(2) = 0;  % at each iteration assume no samplerUpdateOccurred, unless it occurs

        co_proposalFound_samplerUpdateOccurred(1) = 0;  % co_proposalFound = false;
        samplerUpdateIsGreedy = counterAUC < self.SpecDRAM.greedyAdaptationCount.val;

        counterDRS = round(co_AccRate(1));
        if counterDRS > -1, co_proposalFound_samplerUpdateOccurred(1) = 1; end   % co_proposalFound = true

        %*******************************************************************************************************************
        %**************************************** blockProposalAccepted ****************************************************

        % On the very first iteration, this block is (and must be) executed for imageID==1,
        % since it is for the first (starting) point, which is assumed to have been accepted
        % as the first point by the first coarray imageID.

        if co_proposalFound_samplerUpdateOccurred(1) == 1   % blockProposalAccepted: co_proposalAccepted = true
%xxx
            lastState                       = self.Stats.NumFunCall.accepted;
            self.Stats.NumFunCall.accepted  = self.Stats.NumFunCall.accepted + 1;
%xxx
            currentStateWeight              = 0;

            % communicate the accepted logFunc and State from the winning image to master/all images: co_LogFuncState

            % Note: after every adaptive update of the sampler, counterAUP is reset to 1.
            if counterAUP == 0 && samplerUpdateSucceeded
                numFunCallAcceptedLastAdaptation    = numFunCallAcceptedLastAdaptation + 1;
                lastStateWeight                     = 0;
            end

            if self.isFreshRun  % blockFreshDryRun

                self.Chain.ProcessID    (self.Stats.NumFunCall.accepted)    = imageID;
                self.Chain.DelRejStage  (self.Stats.NumFunCall.accepted)    = counterDRS;
                self.Chain.Adaptation   (self.Stats.NumFunCall.accepted)    = chainAdaptationMeasure;
                self.Chain.Weight       (self.Stats.NumFunCall.accepted)    = 0;
                self.Chain.LogFunc      (self.Stats.NumFunCall.accepted)    = co_LogFuncState(1,1);
                self.Chain.State  (1:nd, self.Stats.NumFunCall.accepted)    = co_LogFuncState(2:nd+1,1);

                % find the burnin point

                self.Chain.BurninLoc(self.Stats.NumFunCall.accepted)    = self.getBurninLoc ( self.Stats.NumFunCall.accepted                        ...
                                                                                            , self.Stats.LogFuncMode.val                            ...
                                                                                            , self.Chain.LogFunc(1:self.Stats.NumFunCall.accepted)  ...
                                                                                            ) ;

            else % blockFreshDryRun : in restart mode: determine the correct value of co_proposalFound_samplerUpdateOccurred(1)
                numFunCallAcceptedPlusOne = self.Stats.NumFunCall.accepted + 1;
% disp("numFunCallAcceptedPlusOne = ");
% disp(numFunCallAcceptedPlusOne);
% disp("self.Chain.Count.compact = ");
% disp(self.Chain.Count.compact)
                if numFunCallAcceptedPlusOne == self.Chain.Count.compact
                    self.isFreshRun                                     = true;
                    self.writeOutput();
                    self.isDryRun                                       = ~self.isFreshRun;
                    self.Chain.Weight(numFunCallAcceptedPlusOne)        = 0;
                    SumAccRateSinceStart.acceptedRejected               = self.Chain.MeanAccRate(self.Stats.NumFunCall.accepted) * self.Stats.NumFunCall.acceptedRejected;
                    if delayedRejectionRequested
                        SumAccRateSinceStart.acceptedRejectedDelayed    = self.Chain.MeanAccRate(self.Stats.NumFunCall.accepted) * self.Stats.NumFunCall.acceptedRejectedDelayed;
                    end
                end
                self.Stats.NumFunCall.accepted                          = numFunCallAcceptedPlusOne;
                numFunCallAcceptedPlusOne                               = self.Stats.NumFunCall.accepted + 1;

            end % blockFreshDryRun

            if self.Stats.LogFuncMode.val < self.Chain.LogFunc(self.Stats.NumFunCall.accepted)
                self.Stats.LogFuncMode.val          = self.Chain.LogFunc(self.Stats.NumFunCall.accepted);
                self.Stats.LogFuncMode.Loc.compact  = self.Stats.NumFunCall.accepted;
            end

            SumAccRateSinceStart.acceptedRejected = SumAccRateSinceStart.acceptedRejected + co_AccRate(counterDRS+2);

        else % blockProposalAccepted

            counterDRS                              = self.SpecDRAM.delayedRejectionCount.val;
            SumAccRateSinceStart.acceptedRejected   = SumAccRateSinceStart.acceptedRejected + co_AccRate(counterDRS+2);

        end % blockProposalAccepted

        %**************************************** blockProposalAccepted ****************************************************
        %*******************************************************************************************************************

        counterAUP                              = counterAUP + 1;
        counterPRP                              = counterPRP + 1;
        currentStateWeight                      = currentStateWeight + 1;
        self.Stats.NumFunCall.acceptedRejected  = self.Stats.NumFunCall.acceptedRejected + 1;

        if delayedRejectionRequested
            SumAccRateSinceStart.acceptedRejectedDelayed    = SumAccRateSinceStart.acceptedRejectedDelayed + sum(co_AccRate(2:counterDRS+2));
            self.Stats.NumFunCall.acceptedRejectedDelayed   = self.Stats.NumFunCall.acceptedRejectedDelayed + counterDRS + 1;
        end

        if self.isFreshRun  % these are used for adaptive proposal updating, so they have to be set on every accepted or rejected iteration (excluding delayed rejections)
            self.Chain.MeanAccRate  (self.Stats.NumFunCall.accepted)    = SumAccRateSinceStart.acceptedRejected / self.Stats.NumFunCall.acceptedRejected;
            self.Chain.Weight       (self.Stats.NumFunCall.accepted)    = self.Chain.Weight(self.Stats.NumFunCall.accepted) + 1;
        else
            acceptedRejectedDelayedUnusedRestartMode                    = self.Stats.NumFunCall.acceptedRejectedDelayed;
        end

        if counterPRP == self.SpecBase.progressReportPeriod.val
            counterPRP = 0;
            self.reportProgress();
        end

        %*******************************************************************************************************************
        %******************************************** Proposal Adaptation **************************************************

        if counterAUC < self.SpecDRAM.adaptiveUpdateCount.val && counterAUP  == self.SpecDRAM.adaptiveUpdatePeriod.val    % blockSamplerAdaptation

            co_proposalFound_samplerUpdateOccurred(2) = 1;  % istart = numFunCallAcceptedLastAdaptation % = max( numFunCallAcceptedLastAdaptation , PD%Chain%BurninLoc(PD%Stats%NumFunCall%accepted) ) % this is experimental

            % the order in the following two MUST be preserved as occasionally PD%Stats%NumFunCall%accepted = numFunCallAcceptedLastAdaptation

            dummy = self.Chain.Weight(self.Stats.NumFunCall.accepted);  % needed for the restart mode, not needed in the fresh run
            if self.Stats.NumFunCall.accepted == numFunCallAcceptedLastAdaptation   % no new point has been accepted since last time
                self.Chain.Weight(numFunCallAcceptedLastAdaptation) = currentStateWeight - lastStateWeight;
            else
                self.Chain.Weight(numFunCallAcceptedLastAdaptation) = self.Chain.Weight(numFunCallAcceptedLastAdaptation) - lastStateWeight;
                self.Chain.Weight(self.Stats.NumFunCall.accepted)   = currentStateWeight;   % needed for the restart mode, not needed in the fresh run
            end

            if self.isFreshRun
                meanAccRateSinceStart = self.Chain.MeanAccRate(self.Stats.NumFunCall.accepted);
                if self.SpecBase.restartFileFormat.isBinary
                    fprintf(self.RestartFile.unit, num2str(meanAccRateSinceStart));
                else
                    % self.RestartFile.format
                    fprintf(self.RestartFile.unit, "meanAccRateSinceStart" + "\n" + num2str(meanAccRateSinceStart));
                    self.Proposal.writeRestartFile(); % xxx
                end

            else
                if self.SpecBase.restartFileFormat.isBinary
                % meanAccRateSinceStart = read(self.RestartFile.unit)
                else
                    fgets(self.RestartFile.unit);
                    meanAccRateSinceStart = fgets(self.RestartFile.unit);
                    self.Proposal.readRestartFile();
                end
                SumAccRateSinceStart.acceptedRejected = meanAccRateSinceStart * self.Stats.NumFunCall.acceptedRejected;
            end

            [samplerUpdateSucceeded, adaptationMeasure] = self.Proposal.doAdaptation   ( nd                                                                                        ...
                                                                                            , self.Stats.NumFunCall.accepted - numFunCallAcceptedLastAdaptation + 1                     ...
                                                                                            , self.Chain.State  (1:nd, numFunCallAcceptedLastAdaptation:self.Stats.NumFunCall.accepted) ...
                                                                                            , self.Chain.Weight (numFunCallAcceptedLastAdaptation:self.Stats.NumFunCall.accepted)       ...
                                                                                            , samplerUpdateIsGreedy                                                                     ...
                                                                                            , meanAccRateSinceStart                                                                     ...
                                                                                            , chainAdaptationMeasure                                                                    ...
                                                                                            ) ;

            self.Chain.Weight(self.Stats.NumFunCall.accepted) = dummy;  % needed for the restart mode, not needed in the fresh run
            if self.Stats.NumFunCall.accepted == numFunCallAcceptedLastAdaptation
                %adaptationMeasure = adaptationMeasure + adaptationMeasureDummy;    % this is the worst-case upper-bound
                self.Chain.Adaptation(self.Stats.NumFunCall.accepted)   = self.Chain.Adaptation(self.Stats.NumFunCall.accepted) + adaptationMeasure; % this is the worst-case upper-bound
            else
                %adaptationMeasure = adaptationMeasureDummy
                self.Chain.Adaptation(self.Stats.NumFunCall.accepted)   = adaptationMeasure;
                self.Chain.Weight(numFunCallAcceptedLastAdaptation)     = self.Chain.Weight(numFunCallAcceptedLastAdaptation) + lastStateWeight;
            end
            if samplerUpdateSucceeded
                lastStateWeight = currentStateWeight;   % self.Chain.Weight(self.Stats.NumFunCall.accepted) % informative, do not remove
                numFunCallAcceptedLastAdaptation = self.Stats.NumFunCall.accepted;
            end

            counterAUP = 0;
            counterAUC = counterAUC + 1;
            %if counterAUC == self.SpecDRAM.adaptiveUpdateCount.val, chainAdaptationMeasure = 0; end

        end % blockSamplerAdaptation

        %******************************************** Proposal Adaptation **************************************************
        %*******************************************************************************************************************


        %*******************************************************************************************************************
        %************************************************* output write ****************************************************

        % if new point has been sampled, write the previous sampled point to output file

        if self.isFreshRun && co_proposalFound_samplerUpdateOccurred(1) == 1 && lastState > 0   % blockOutputWrite

            if self.SpecBase.chainFileFormat.isCompact
                fprintf(self.ChainFile.unit , self.ChainFile.format                     ...
                                            , self.Chain.ProcessID  (lastState)         ...
                                            , self.Chain.DelRejStage(lastState)         ...
                                            , self.Chain.MeanAccRate(lastState)         ...
                                            , self.Chain.Adaptation (lastState)         ...
                                            , self.Chain.BurninLoc  (lastState)         ...
                                            , self.Chain.Weight     (lastState)         ...
                                            , self.Chain.LogFunc    (lastState)         ...
                                            , self.Chain.State (1:nd,lastState)         ...
                                            ) ;
            elseif self.SpecBase.chainFileFormat.isBinary

                fprintf(self.ChainFile.unit , self.ChainFile.format                     ...
                                            , self.Chain.ProcessID  (lastState)         ...
                                            , self.Chain.DelRejStage(lastState)         ...
                                            , self.Chain.MeanAccRate(lastState)         ...
                                            , self.Chain.Adaptation (lastState)         ...
                                            , self.Chain.BurninLoc  (lastState)         ...
                                            , 1                                         ...
                                            , self.Chain.LogFunc    (lastState)         ...
                                            , self.Chain.State (1:nd,lastState)         ...
                                            ) ;
            elseif self.SpecBase.chainFileFormat.isVerbose

                for j = 1 : self.Chain.Weight(lastState)
                fprintf(self.ChainFile.unit , self.ChainFile.format                     ...
                                            , self.Chain.ProcessID  (lastState)         ...
                                            , self.Chain.DelRejStage(lastState)         ...
                                            , self.Chain.MeanAccRate(lastState)         ...
                                            , self.Chain.Adaptation (lastState)         ...
                                            , self.Chain.BurninLoc  (lastState)         ...
                                            , 1                                         ...
                                            , self.Chain.LogFunc    (lastState)         ...
                                            , self.Chain.State (1:nd,lastState)         ...
                                            ) ;
                end

            end

        end % blockOutputWrite

        %************************************************* output write ****************************************************
        %*******************************************************************************************************************

        %*******************************************************************************************************************
        %********************************************* last output write ***************************************************

        if self.Stats.NumFunCall.accepted == self.SpecMCMC.chainSize.val    % blockLastSample: co_missionAccomplished = true

            % on 3 images Windows, sustituting co_missionAccomplished with the following leads to 10% less communication overhead for 1D Gaussian example
            % on 3 images Linux  , sustituting co_missionAccomplished with the following leads to 16% less communication overhead for 1D Gaussian example
            % on 5 images Linux  , sustituting co_missionAccomplished with the following leads to 11% less communication overhead for 1D Gaussian example

            co_proposalFound_samplerUpdateOccurred(1) = -1; % equivalent to co_missionAccomplished = true
            inverseProgressReportPeriod = 1 / (self.Stats.NumFunCall.acceptedRejected - numFunCallAcceptedRejectedLastReport);

            if self.isFreshRun  % blockLastOutputWrite

                % write the last sampled point to the output files

                lastState = self.Stats.NumFunCall.accepted;

                if self.SpecBase.chainFileFormat.isCompact || self.SpecBase.chainFileFormat.isVerbose
                    fprintf(self.ChainFile.unit , self.ChainFile.format                 ...
                                                , self.Chain.ProcessID  (lastState)     ...
                                                , self.Chain.DelRejStage(lastState)     ...
                                                , self.Chain.MeanAccRate(lastState)     ...
                                                , self.Chain.Adaptation (lastState)     ...
                                                , self.Chain.BurninLoc  (lastState)     ...
                                                , 1                                     ...
                                                , self.Chain.LogFunc    (lastState)     ...
                                                , self.Chain.State (1:nd,lastState)     ...
                                                ) ;
                elseif self.SpecBase.ChainFileFormat.isBinary
                    fprintf(self.ChainFile.unit , self.ChainFile.format                 ...
                                                , self.Chain.ProcessID  (lastState)     ...
                                                , self.Chain.DelRejStage(lastState)     ...
                                                , self.Chain.MeanAccRate(lastState)     ...
                                                , self.Chain.Adaptation (lastState)     ...
                                                , self.Chain.BurninLoc  (lastState)     ...
                                                , 1                                     ...
                                                , self.Chain.LogFunc    (lastState)     ...
                                                , self.Chain.State (1:nd,lastState)     ...
                                                ) ;
                end

            end % blockLastOutputWrite

            self.reportProgress();

        end % blockLastSample

        %********************************************* last output write ***************************************************
        %*******************************************************************************************************************

        if co_proposalFound_samplerUpdateOccurred(1) == -1, break; end % loopMarkovChain: we are done: co_missionAccomplished = true

        co_AccRate(1) = -1; % counterDRS at which new proposal is accepted
        maxLogFuncRejectedProposal = Constants.NEGINF_RK;

        for counterDRS = 0 : self.SpecDRAM.delayedRejectionCount.val    % loopDelayedRejection

            co_LogFuncState(2:nd+1,counterDRS+2)    = self.Proposal.getNew  ( nd                                    ...
                                                                            , counterDRS                            ...
                                                                            , co_LogFuncState(2:nd+1, counterDRS+1) ...
                                                                            ) ;

            uniformRnd = rand;  % only for the purpose of restart mode reproducibility

            if self.isFreshRun || numFunCallAcceptedPlusOne == self.Chain.Count.compact

                self.Timer.toc();
                co_LogFuncState(1, counterDRS+2)    = getLogFunc(co_LogFuncState(2:nd+1, counterDRS+2));
                self.Timer.toc();
                self.Stats.avgTimePerFunCalInSec    = self.Stats.avgTimePerFunCalInSec + self.Timer.delta;

                % accept or reject the proposed state

                if co_LogFuncState(1,counterDRS+2)  >= co_LogFuncState(1,1) % accept the proposed state

                    co_AccRate(counterDRS+2)    = 1.0;
                    co_AccRate(1)               = counterDRS;
                    co_LogFuncState(1:nd+1,1)   = co_LogFuncState(1:nd+1,counterDRS+2);
                    break   % loopDelayedRejection

                elseif co_LogFuncState(1,counterDRS+2) < maxLogFuncRejectedProposal % reject the proposed state. This step should be reachable only when delayedRejectionCount > 0

                    co_AccRate(counterDRS+2) = 0.0;   % proposal rejected

                else % accept with probability co_AccRate

                    if counterDRS == 0  % This should be equivalent to maxLogFuncRejectedProposal == NEGINF_RK
                        co_AccRate(counterDRS+2)    = exp( co_LogFuncState(1,counterDRS+2) - co_LogFuncState(1,1) );
                    else % ensure no arithmetic overflow/underflow. ATTN: co_LogFuncState(0,-1) > co_LogFuncState(0,counterDRS) > maxLogFuncRejectedProposal
                        co_AccRate(counterDRS+2)    = exp( getLogSubExp(co_LogFuncState(1,counterDRS+2) , maxLogFuncRejectedProposal)       ...
                                                         - getLogSubExp(co_LogFuncState(1,1)            , maxLogFuncRejectedProposal) );
                    end
                    if uniformRnd < co_AccRate(counterDRS+2)    % accept the proposed state
                        co_AccRate(1)               = counterDRS;
                        co_LogFuncState(1:nd+1,1)   = co_LogFuncState(1:nd+1,counterDRS+2);
                        break   % loopDelayedRejection
                    end

                end

                maxLogFuncRejectedProposal = max( maxLogFuncRejectedProposal, co_LogFuncState(1,counterDRS+2) );

            else

                if currentStateWeight == self.Chain.Weight(self.Stats.NumFunCall.accepted) && counterDRS == self.Chain.DelRejStage(numFunCallAcceptedPlusOne)
                    co_AccRate(1)               = counterDRS;
                    co_LogFuncState(   1, 1)    = self.Chain.LogFunc(numFunCallAcceptedPlusOne);
                    co_LogFuncState(2:nd+1,1)   = self.Chain.State  (1:nd,numFunCallAcceptedPlusOne);
                    break   % loopDelayedRejection
                end

            end

        end % loopDelayedRejection

        %***************************************** end Common block between all images *************************************
        %*******************************************************************************************************************

    end % loopMarkovChain

    self.Timer.toc();

    %********************************************** end of loopMarkovChain *************************************************
    %***********************************************************************************************************************
    %***********************************************************************************************************************

    self.Chain.Count.resized = self.Stats.NumFunCall.accepted;
    self.Chain.Count.compact = self.Stats.NumFunCall.accepted;
    self.Chain.Count.verbose = self.Stats.NumFunCall.acceptedRejected;

    if noDelayedRejectionRequested
        self.Stats.NumFunCall.acceptedRejectedDelayed   = self.Stats.NumFunCall.acceptedRejected;
        SumAccRateSinceStart.acceptedRejectedDelayed    = SumAccRateSinceStart.acceptedRejected;
    end

    if self.Image.isMaster

        if self.SpecDRAM.burninAdaptationMeasure.val > 0.999999
            self.Stats.AdaptationBurninLoc.compact = 1;
            self.Stats.AdaptationBurninLoc.verbose = 1;
        else
            self.Stats.AdaptationBurninLoc.compact = self.Stats.NumFunCall.accepted;
            self.Stats.AdaptationBurninLoc.verbose = self.Stats.NumFunCall.acceptedRejected - self.Chain.Weight(self.Stats.NumFunCall.accepted) + 1;
            for i = self.Stats.NumFunCall.accepted-1 : -1 : 1
                if self.Chain.Adaptation(i) > self.SpecDRAM.burninAdaptationMeasure.val, break; end
                self.Stats.AdaptationBurninLoc.compact = self.Stats.AdaptationBurninLoc.compact - 1;
                self.Stats.AdaptationBurninLoc.verbose = self.Stats.AdaptationBurninLoc.verbose - self.Chain.Weight(i);
            end
        end

        self.Stats.BurninLoc.compact          = self.Chain.BurninLoc(self.Stats.NumFunCall.accepted);
        self.Stats.BurninLoc.verbose          = sum(self.Chain.Weight(1:self.Stats.BurninLoc.compact-1)) + 1;
        self.Stats.LogFuncMode.Crd            = self.Chain.State(1:nd,self.Stats.LogFuncMode.Loc.compact);
        self.Stats.LogFuncMode.Loc.verbose    = sum(self.Chain.Weight(1:self.Stats.LogFuncMode.Loc.compact-1)) + 1;

        fclose(self.RestartFile.unit);
        fclose(self.ChainFile.unit);
    end

    self.Stats.avgCommTimePerFunCall                    = 0;
    self.Stats.NumFunCall.acceptedRejectedDelayedUnused = self.Stats.NumFunCall.acceptedRejectedDelayed;
    self.Stats.avgTimePerFunCalInSec                    = self.Stats.avgTimePerFunCalInSec / (self.Stats.NumFunCall.acceptedRejectedDelayedUnused - acceptedRejectedDelayedUnusedRestartMode);

    %***********************************************************************************************************************
    %***********************************************************************************************************************

end
