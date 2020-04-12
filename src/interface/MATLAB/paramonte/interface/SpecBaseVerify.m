%**********************************************************************************************************************************
%**********************************************************************************************************************************
%
%  ParaMonte: plain powerful parallel Monte Carlo library.
%
%  Copyright (C) 2012-present, The Computational Data Science Lab
%
%  This file is part of ParaMonte library. 
%
%  ParaMonte is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, version 3 of the License.
%
%  ParaMonte is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
%
%**********************************************************************************************************************************
%**********************************************************************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SpecBase specification type-checking class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef SpecBaseVerify

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        delim = ",";
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %function self = SpecBaseVerify()
        %end

        function result = sampleSize(self,sampleSize)
            result = verifySpec(sampleSize,self.objectName+"sampleSize","real",delim);
        end

        function result = randomSeed(self,randomSeed)
            result = verifySpec(randomSeed,self.objectName+"randomSeed","integer",delim);
        end

        function result = description(self,description)
            result = verifySpec(description,self.objectName+"description","string",delim);
        end

        function result = outputFileName(self,outputFileName)
            result = verifySpec(outputFileName,self.objectName+"outputFileName","string",delim);
        end

        function result = outputDelimiter(self,outputDelimiter)
            result = verifySpec(outputDelimiter,self.objectName+"outputDelimiter","string",delim);
        end

        function result = chainFileFormat(self,chainFileFormat)
            result = verifySpec(chainFileFormat,self.objectName+"chainFileFormat","string",delim);
        end

        function result = variableNameList(self,variableNameList)
            result = verifySpec(variableNameList,self.objectName+"variableNameList","string",delim);
        end

        function result = restartFileFormat(self,restartFileFormat)
            result = verifySpec(restartFileFormat,self.objectName+"restartFileFormat","string",delim);
        end

        function result = outputColumnWidth(self,outputColumnWidth)
            result = verifySpec(outputColumnWidth,self.objectName+"outputColumnWidth","integer",delim);
        end

        function result = outputRealPrecision(self,outputRealPrecision)
            result = verifySpec(outputRealPrecision,self.objectName+"outputRealPrecision","integer",delim);
        end

        function result = silentModeRequested(self,silentModeRequested)
            result = verifySpec(silentModeRequested,self.objectName+"silentModeRequested","logical",delim);
        end

        function result = domainLowerLimitVec(self,domainLowerLimitVec)
            result = verifySpec(domainLowerLimitVec,self.objectName+"domainLowerLimitVec","real",delim);
        end

        function result = domainUpperLimitVec(self,domainUpperLimitVec)
            result = verifySpec(domainUpperLimitVec,self.objectName+"domainUpperLimitVec","real",delim);
        end

        function result = parallelizationModel(self,parallelizationModel)
            result = verifySpec(parallelizationModel,self.objectName+"parallelizationModel","string",delim);
        end

        function result = progressReportPeriod(self,progressReportPeriod)
            result = verifySpec(progressReportPeriod,self.objectName+"progressReportPeriod","integer",delim);
        end

        function result = targetAcceptanceRate(self,targetAcceptanceRate)
            result = verifySpec(targetAcceptanceRate,self.objectName+"targetAcceptanceRate","real",delim);
        end

        function result = mpiFinalizeRequested(self,mpiFinalizeRequested)
            result = verifySpec(mpiFinalizeRequested,self.objectName+"mpiFinalizeRequested","logical",delim);
        end

        function result = maxNumDomainCheckToWarn(self,maxNumDomainCheckToWarn)
            result = verifySpec(maxNumDomainCheckToWarn,self.objectName+"maxNumDomainCheckToWarn","integer",delim);
        end

        function result = maxNumDomainCheckToStop(self,maxNumDomainCheckToStop)
            result = verifySpec(maxNumDomainCheckToStop,self.objectName+"maxNumDomainCheckToStop","integer",delim);
        end

        function result = interfaceType(self)
            result = "interfaceType=" + "'MATLAB " + version + "'" + self.delim
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecBase