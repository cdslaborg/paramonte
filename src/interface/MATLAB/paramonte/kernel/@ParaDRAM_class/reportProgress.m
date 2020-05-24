function reportProgress(self)

    global  timeElapsedUntilLastReportInSeconds SumAccRateSinceStart sumAccRateLastReport inverseProgressReportPeriod    ...
            numFunCallAcceptedRejectedLastReport numFunCallAccepted_dummy

    if self.isFreshRun
        self.Timer.toc();
        timeElapsedSinceLastReportInSeconds         = self.Timer.total - timeElapsedUntilLastReportInSeconds;
        timeElapsedUntilLastReportInSeconds         = self.Timer.total;
        meanAccRateSinceStart               = SumAccRateSinceStart.acceptedRejected / self.Stats.NumFunCall.acceptedRejected;
        meanAccRateSinceLastReport          = (SumAccRateSinceStart.acceptedRejected - sumAccRateLastReport) * inverseProgressReportPeriod;
        estimatedTimeToFinishInSeconds      = (self.SpecMCMC.chainSize.val - self.Stats.NumFunCall.accepted) * timeElapsedUntilLastReportInSeconds / self.Stats.NumFunCall.accepted;
       
        % self.TimeFile.format
        fprintf ( self.TimeFile.unit    , self.TimeFile.format                      ...
                                        , self.Stats.NumFunCall.acceptedRejected    ...
                                        , self.Stats.NumFunCall.accepted            ...
                                        , meanAccRateSinceStart                     ...
                                        , meanAccRateSinceLastReport                ...
                                        , timeElapsedSinceLastReportInSeconds       ...
                                        , timeElapsedUntilLastReportInSeconds       ...
                                        , estimatedTimeToFinishInSeconds            ...
                                        ) ;
    else

        Record.value    = fgets(self.TimeFile.unit);
        Record.Parts    = split(strtrim(Record.value), self.SpecBase.outputDelimiter.val);
        
        numFunCallAcceptedRejectedLastReport    = str2num(Record.Parts{1});
        numFunCallAccepted_dummy                = str2num(Record.Parts{2});
        meanAccRateSinceStart                   = str2num(Record.Parts{3});
        meanAccRateSinceLastReport              = str2num(Record.Parts{4});
        timeElapsedSinceLastReportInSeconds     = str2num(Record.Parts{5});
        timeElapsedUntilLastReportInSeconds     = str2num(Record.Parts{6});
        estimatedTimeToFinishInSeconds          = str2num(Record.Parts{7});

        SumAccRateSinceStart.acceptedRejected   = meanAccRateSinceStart * numFunCallAcceptedRejectedLastReport;

    end

    % report progress in the standard output

    formatStr   = "            %14d / %8d" + "   " + "%16.4f / %6.4f" + "   " + "%14.4f / %8.4f" ;
    if self.Image.isFirst
        fprintf ( 1, repmat('\b', 1, 93));
        fprintf ( 1, formatStr  , self.Stats.NumFunCall.accepted                                                    ...
                                , self.Stats.NumFunCall.acceptedRejected                                            ...
                                , meanAccRateSinceLastReport                                                        ...
                                , SumAccRateSinceStart.acceptedRejected / self.Stats.NumFunCall.acceptedRejected    ...
                                , timeElapsedUntilLastReportInSeconds                                               ...
                                , estimatedTimeToFinishInSeconds                                                    ...
                                ) ;
    end

    numFunCallAcceptedRejectedLastReport    = self.Stats.NumFunCall.acceptedRejected;
    sumAccRateLastReport                    = SumAccRateSinceStart.acceptedRejected;

end
