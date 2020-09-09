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

function reportProgress(self)

%    global  timeElapsedUntilLastReportInSeconds SumAccRateSinceStart sumAccRateLastReport inverseProgressReportPeriod    ...
%            numFunCallAcceptedRejectedLastReport numFunCallAccepted_dummy

    if self.isFreshRun

        self.Timer.toc();
        timeElapsedSinceLastReportInSeconds                         = self.Timer.total - self.Proposal.Global.timeElapsedUntilLastReportInSeconds;
        self.Proposal.Global.timeElapsedUntilLastReportInSeconds    = self.Timer.total;
        meanAccRateSinceStart                                       = self.Proposal.Global.SumAccRateSinceStart.acceptedRejected / self.Stats.NumFunCall.acceptedRejected;
        meanAccRateSinceLastReport                                  = (self.Proposal.Global.SumAccRateSinceStart.acceptedRejected - self.Proposal.Global.sumAccRateLastReport) * self.Proposal.Global.inverseProgressReportPeriod;
        estimatedTimeToFinishInSeconds                              = (self.SpecMCMC.chainSize.val - self.Stats.NumFunCall.accepted) * self.Proposal.Global.timeElapsedUntilLastReportInSeconds / self.Stats.NumFunCall.accepted;
       
        fprintf ( self.TimeFile.unit    , self.TimeFile.format                                      ...
                                        , self.Stats.NumFunCall.acceptedRejected                    ...
                                        , self.Stats.NumFunCall.accepted                            ...
                                        , meanAccRateSinceStart                                     ...
                                        , meanAccRateSinceLastReport                                ...
                                        , timeElapsedSinceLastReportInSeconds                       ...
                                        , self.Proposal.Global.timeElapsedUntilLastReportInSeconds  ...
                                        , estimatedTimeToFinishInSeconds                            ...
                                        ) ;

    else

        Record.value = fgets(self.TimeFile.unit);
        if Record.value ~= -1
            Record.Parts    = split(strtrim(Record.value), self.SpecBase.outputDelimiter.val);
            self.Proposal.Global.numFunCallAcceptedRejectedLastReport   = str2num(Record.Parts{1});
            self.Proposal.Global.numFunCallAccepted_dummy               = str2num(Record.Parts{2});
            meanAccRateSinceStart                                       = str2num(Record.Parts{3});
            meanAccRateSinceLastReport                                  = str2num(Record.Parts{4});
            timeElapsedSinceLastReportInSeconds                         = str2num(Record.Parts{5});
            self.Proposal.Global.timeElapsedUntilLastReportInSeconds    = str2num(Record.Parts{6});
            estimatedTimeToFinishInSeconds                              = str2num(Record.Parts{7});
            self.Proposal.Global.SumAccRateSinceStart.acceptedRejected  = meanAccRateSinceStart * self.Proposal.Global.numFunCallAcceptedRejectedLastReport;
        end

    end

    % report progress in the standard output

    formatStr   = "            %14d / %8d" + "   " + "%16.4f / %6.4f" + "   " + "%11.2f / %11.2f" ;
    if self.Image.isFirst
        fprintf ( 1, repmat('\b', 1, 93));
        fprintf ( 1, formatStr  , self.Stats.NumFunCall.accepted                                                                        ...
                                , self.Stats.NumFunCall.acceptedRejected                                                                ...
                                , meanAccRateSinceLastReport                                                                            ...
                                , self.Proposal.Global.SumAccRateSinceStart.acceptedRejected / self.Stats.NumFunCall.acceptedRejected   ...
                                , self.Proposal.Global.timeElapsedUntilLastReportInSeconds                                              ...
                                , estimatedTimeToFinishInSeconds                                                                        ...
                                ) ;
    end

    self.Proposal.Global.numFunCallAcceptedRejectedLastReport   = self.Stats.NumFunCall.acceptedRejected;
    self.Proposal.Global.sumAccRateLastReport                   = self.Proposal.Global.SumAccRateSinceStart.acceptedRejected;

end
