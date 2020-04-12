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

classdef SpecBaseVerify < SpecVerification

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        function self = SpecBaseVerify(objectName)
            self.objectName = objectName;
        end

        function result = sampleSize(self,sampleSize)
            result = self.verifySpec(sampleSize,"real");
        end

        function result = randomSeed(self,randomSeed)
            result = self.verifySpec(randomSeed,"integer");
        end

        function result = description(self,description)
            result = self.verifySpec(description,"string");
        end

        function result = outputFileName(self,outputFileName)
            result = self.verifySpec(outputFileName,"string");
        end

        function result = outputDelimiter(self,outputDelimiter)
            result = self.verifySpec(outputDelimiter,"string");
        end

        function result = chainFileFormat(self,chainFileFormat)
            result = self.verifySpec(chainFileFormat,"string");
        end

        function result = variableNameList(self,variableNameList)
            result = self.verifySpec(variableNameList,"string");
        end

        function result = restartFileFormat(self,restartFileFormat)
            result = self.verifySpec(restartFileFormat,"string");
        end

        function result = outputColumnWidth(self,outputColumnWidth)
            result = self.verifySpec(outputColumnWidth,"integer");
        end

        function result = outputRealPrecision(self,outputRealPrecision)
            result = self.verifySpec(outputRealPrecision,"integer");
        end

        function result = silentModeRequested(self,silentModeRequested)
            result = self.verifySpec(silentModeRequested,"logical");
        end

        function result = domainLowerLimitVec(self,domainLowerLimitVec)
            result = self.verifySpec(domainLowerLimitVec,"real");
        end

        function result = domainUpperLimitVec(self,domainUpperLimitVec)
            result = self.verifySpec(domainUpperLimitVec,"real");
        end

        function result = parallelizationModel(self,parallelizationModel)
            result = self.verifySpec(parallelizationModel,"string");
        end

        function result = progressReportPeriod(self,progressReportPeriod)
            result = self.verifySpec(progressReportPeriod,"integer");
        end

        function result = targetAcceptanceRate(self,targetAcceptanceRate)
            result = self.verifySpec(targetAcceptanceRate,"real");
        end

        function result = mpiFinalizeRequested(self,mpiFinalizeRequested)
            result = self.verifySpec(mpiFinalizeRequested,"logical");
        end

        function result = maxNumDomainCheckToWarn(self,maxNumDomainCheckToWarn)
            result = self.verifySpec(maxNumDomainCheckToWarn,"integer");
        end

        function result = maxNumDomainCheckToStop(self,maxNumDomainCheckToStop)
            result = self.verifySpec(maxNumDomainCheckToStop,"integer");
        end

        function result = interfaceType(self)
            result = "interfaceType=" + "'MATLAB " + version + "'" + self.delim;
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecBase