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
%%%% SpecDRAM specification type-checking class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef SpecDRAM

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        delim = ",";
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %function self = SpecDRAM()
        %end

        function result = scaleFactor(self,scaleFactor)
            result = verifySpec(scaleFactor,self.objectName+"scaleFactor","string",delim);
        end

        function result = proposalModel(self,proposalModel)
            result = verifySpec(proposalModel,self.objectName+"proposalModel","string",delim);
        end

        function result = proposalStartCovMat(self,proposalStartCovMat)
            result = verifySpec(proposalStartCovMat,self.objectName+"proposalStartCovMat","real",delim);
        end

        function result = proposalStartCorMat(self,proposalStartCorMat)
            result = verifySpec(proposalStartCorMat,self.objectName+"proposalStartCorMat","real",delim);
        end

        function result = proposalStartStdVec(self,proposalStartStdVec)
            result = verifySpec(proposalStartStdVec,self.objectName+"proposalStartStdVec","real",delim);
        end

        function result = adaptiveUpdateCount(self,adaptiveUpdateCount)
            result = verifySpec(adaptiveUpdateCount,self.objectName+"adaptiveUpdateCount","integer",delim);
        end

        function result = adaptiveUpdatePeriod(self,adaptiveUpdatePeriod)
            result = verifySpec(adaptiveUpdatePeriod,self.objectName+"adaptiveUpdatePeriod","integer",delim);
        end

        function result = greedyAdaptationCount(self,greedyAdaptationCount)
            result = verifySpec(greedyAdaptationCount,self.objectName+"greedyAdaptationCount","integer",delim);
        end

        function result = delayedRejectionCount(self,delayedRejectionCount)
            result = verifySpec(delayedRejectionCount,self.objectName+"delayedRejectionCount","integer",delim);
        end

        function result = burninAdaptationMeasure(self,burninAdaptationMeasure)
            result = verifySpec(burninAdaptationMeasure,self.objectName+"burninAdaptationMeasure","integer",delim);
        end

        function result = delayedRejectionScaleFactorVec(self,delayedRejectionScaleFactorVec)
            result = verifySpec(delayedRejectionScaleFactorVec,self.objectName+"delayedRejectionScaleFactorVec","real",delim);
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecDRAM