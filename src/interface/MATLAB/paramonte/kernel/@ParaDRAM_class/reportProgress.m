%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it 
%   under the terms of the GNU Lesser General Public License as published 
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see, 
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms, 
%   if you use any parts of this library for any purposes, 
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
        if Record.value ~= -1
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

    end

    % report progress in the standard output

    formatStr   = "            %14d / %8d" + "   " + "%16.4f / %6.4f" + "   " + "%11.2f / %11.2f" ;
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
