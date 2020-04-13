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
%%%% SpecMCMC specification type-checking class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SpecMCMCVerify < SpecVerification

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        function self = SpecMCMCVerify(objectName)
            self.objectName = objectName;
        end

        function result = chainSize(self,chainSize)
            result = self.verifySpec(chainSize,"integer");
        end

        function result = startPointVec(self,startPointVec)
            result = self.verifySpec(startPointVec,"real");
        end

        function result = sampleRefinementCount(self,sampleRefinementCount)
            result = self.verifySpec(sampleRefinementCount,"integer");
        end

        function result = sampleRefinementMethod(self,sampleRefinementMethod)
            result = self.verifySpec(sampleRefinementMethod,"string");
        end

        function result = randomStartPointRequested(self,randomStartPointRequested)
            result = self.verifySpec(randomStartPointRequested,"logical");
        end

        function result = randomStartPointDomainLowerLimitVec(self,randomStartPointDomainLowerLimitVec)
            result = self.verifySpec(randomStartPointDomainLowerLimitVec,"real");
        end

        function result = randomStartPointDomainUpperLimitVec(self,randomStartPointDomainUpperLimitVec)
            result = self.verifySpec(randomStartPointDomainUpperLimitVec,"real");
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecDRAMVerify