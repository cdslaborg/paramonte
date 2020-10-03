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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SpecBase specification type-checking class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SpecBaseVerify < SpecVerification

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        function self = SpecBaseVerify(objectName,ndim)
            self.objectName = objectName;
            self.ndim = ndim;
        end

        function result = sampleSize(self,sampleSize)
            result = self.verifySpec(sampleSize,"real",1);
        end

        function result = randomSeed(self,randomSeed)
            result = self.verifySpec(randomSeed,"integer",1);
        end

        function result = description(self,description)
            result = self.verifySpec(description,"string",1);
        end

        function result = outputFileName(self,outputFileName)
            result = self.verifySpec(outputFileName,"string",1);
        end

        function result = outputDelimiter(self,outputDelimiter)
            result = self.verifySpec(outputDelimiter,"string",1);
        end

        function result = chainFileFormat(self,chainFileFormat)
            result = self.verifySpec(chainFileFormat,"string",1);
        end

        function result = variableNameList(self,variableNameList)
            result = self.verifySpec(variableNameList,"string",self.ndim);
        end

        function result = restartFileFormat(self,restartFileFormat)
            result = self.verifySpec(restartFileFormat,"string",1);
        end

        function result = outputColumnWidth(self,outputColumnWidth)
            result = self.verifySpec(outputColumnWidth,"integer",1);
        end

        function result = overwriteRequested(self,overwriteRequested)
            result = self.verifySpec(overwriteRequested,"logical",1);
        end

        function result = outputRealPrecision(self,outputRealPrecision)
            result = self.verifySpec(outputRealPrecision,"integer",1);
        end

        function result = silentModeRequested(self,silentModeRequested)
            result = self.verifySpec(silentModeRequested,"logical",1);
        end

        function result = domainLowerLimitVec(self,domainLowerLimitVec)
            result = self.verifySpec(domainLowerLimitVec,"real",self.ndim);
        end

        function result = domainUpperLimitVec(self,domainUpperLimitVec)
            result = self.verifySpec(domainUpperLimitVec,"real",self.ndim);
        end

        function result = parallelizationModel(self,parallelizationModel)
            result = self.verifySpec(parallelizationModel,"string",1);
        end

        function result = progressReportPeriod(self,progressReportPeriod)
            result = self.verifySpec(progressReportPeriod,"integer",1);
        end

        function result = targetAcceptanceRate(self,targetAcceptanceRate)
            result = self.verifySpec(targetAcceptanceRate,"real",2);
        end

        function result = mpiFinalizeRequested(self,mpiFinalizeRequested)
            result = self.verifySpec(mpiFinalizeRequested,"logical",1);
        end

        function result = maxNumDomainCheckToWarn(self,maxNumDomainCheckToWarn)
            result = self.verifySpec(maxNumDomainCheckToWarn,"integer",1);
        end

        function result = maxNumDomainCheckToStop(self,maxNumDomainCheckToStop)
            result = self.verifySpec(maxNumDomainCheckToStop,"integer",1);
        end

        function result = interfaceType(self)
            result = "interfaceType=" + "'MATLAB " + version + "'" + self.delim;
        end

        function result = systemInfoFilePath(self,systemInfoFilePath)
            if isempty(systemInfoFilePath)
                result = "";
            else
                result = "systemInfoFilePath='" + systemInfoFilePath + "'" + self.delim;
            end
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end