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

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecBase