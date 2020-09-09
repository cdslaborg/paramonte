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

function iac = getBatchMeansIAC(np,Point,Weight,batchSize)

    % note that np must be large enough to get a meaningful answer
    
    FUNCTION_NAME   = "getBatchMeansIAC()";
    
    CumSumWeight    = zeros(1, np);
    
    if ~isempty(Weight)
        CumSumWeight        = cumsum(Weight);
    else
        CumSumWeight(np)    = np;
    end

    % compute batch size and count
    
    if ~isempty(batchSize)
        sampleSize = batchSize;
    else
        sampleSize = floor( CumSumWeight(np)^(0.666666666666666) );
    end
    sampleSizeInverse   = 1 / sampleSize;
    sampleCount         = floor( CumSumWeight(np) / sampleSize );   % xxx: round()
    npEffective         = sampleCount * sampleSize;

    if sampleCount < 2
        iac = 1;
        return
    end

    BatchMean = zeros(1, sampleCount);

    % improvement: iterate from the end to the beginning of the chain to ignore initial points instead of the last points.
    % this would be beneficial for MCMC samples

    % compute the Batch-Means avergage and variance, also average of Point
    
    avgPoint    = 0;
    ipStart     = 1;
    ipEnd       = 0;
    if ~isempty(Weight)
        ip                  = 1;
        isample             = 1;
        ipVerbose           = 0;
        currentSampleEndLoc = sampleSize;
        BatchMean(isample)  = 0;
        while true  % loopOverWeight
            ipVerbose = ipVerbose + 1;
            if ipVerbose > CumSumWeight(ip), ip = ip + 1; end
            if ipVerbose > currentSampleEndLoc  % we are done with the current batch
                avgPoint            = avgPoint + BatchMean(isample);
                BatchMean(isample)  = BatchMean(isample) * sampleSizeInverse;
                if ipVerbose > npEffective, break; end  % exit loopOverWeight: condition equivalent to currentSampleEndLoc == npEffective
                currentSampleEndLoc = currentSampleEndLoc + sampleSize;
                isample             = isample + 1;
                BatchMean(isample)  = 0;
            end

            BatchMean(isample)  = BatchMean(isample) + Point(ip);
        end % loopOverWeight
    else    % there is no weight
        for isample = 1 : sampleCount
            BatchMean(isample)  = 0;
            ipEnd               = ipEnd + sampleSize;
            for ip = ipStart : ipEnd
                BatchMean(isample)  = BatchMean(isample) + Point(ip);
            end
            ipStart             = ipEnd + 1;
            avgPoint            = avgPoint + BatchMean(isample);
            BatchMean(isample)  = BatchMean(isample) * sampleSizeInverse;
        end
    end
    avgBatchMean = sum(BatchMean) / sampleCount;
    varBatchMean = sum((BatchMean - avgBatchMean).^2) / (sampleCount-1);
    avgPoint = avgPoint / npEffective;

    % compute the variance of Point
    
    varPoint = 0;
    if ~isempty(Weight)
        ip          = 1;
        ipVerbose   = 0;
        diffSquared = (Point(ip) - avgPoint)^2;
        while true% loopComputeVarPoint
            ipVerbose = ipVerbose + 1;
            if ipVerbose > npEffective, break; end% exit loopComputeVarPoint
            if ipVerbose > CumSumWeight(ip)
                ip          = ip + 1 ; % by definition, ip never become > np, otherwise it's a disastrous coding bug
                diffSquared = (Point(ip) - avgPoint)^2;
            end
            varPoint = varPoint + diffSquared ;
        end % loopComputeVarPoint
    else
        for ip = 1 : npEffective
            varPoint = varPoint + (Point(ip) - avgPoint)^2;
        end
    end
    varPoint = varPoint / (npEffective-1);

    % compute the IAC

    iac = sampleSize * varBatchMean / varPoint;

end