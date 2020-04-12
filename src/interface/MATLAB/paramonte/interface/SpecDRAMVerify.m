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

classdef SpecDRAMVerify < SpecVerification

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        function self = SpecDRAMVerify(objectName)
            self.objectName = objectName;
        end

        function result = scaleFactor(self,scaleFactor)
            result = self.verifySpec(scaleFactor,"string");
        end

        function result = proposalModel(self,proposalModel)
            result = self.verifySpec(proposalModel,"string");
        end

        function result = proposalStartCovMat(self,proposalStartCovMat)
            result = self.verifySpec(proposalStartCovMat,"real");
        end

        function result = proposalStartCorMat(self,proposalStartCorMat)
            result = self.verifySpec(proposalStartCorMat,"real");
        end

        function result = proposalStartStdVec(self,proposalStartStdVec)
            result = self.verifySpec(proposalStartStdVec,"real");
        end

        function result = adaptiveUpdateCount(self,adaptiveUpdateCount)
            result = self.verifySpec(adaptiveUpdateCount,"integer");
        end

        function result = adaptiveUpdatePeriod(self,adaptiveUpdatePeriod)
            result = self.verifySpec(adaptiveUpdatePeriod,"integer");
        end

        function result = greedyAdaptationCount(self,greedyAdaptationCount)
            result = self.verifySpec(greedyAdaptationCount,"integer");
        end

        function result = delayedRejectionCount(self,delayedRejectionCount)
            result = self.verifySpec(delayedRejectionCount,"integer");
        end

        function result = burninAdaptationMeasure(self,burninAdaptationMeasure)
            result = self.verifySpec(burninAdaptationMeasure,"integer");
        end

        function result = delayedRejectionScaleFactorVec(self,delayedRejectionScaleFactorVec)
            result = self.verifySpec(delayedRejectionScaleFactorVec,"real");
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecDRAMVerify