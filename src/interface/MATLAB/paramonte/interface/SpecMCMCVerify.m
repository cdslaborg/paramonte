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


classdef SpecMCMC

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        delim = ",";
        objectName = "ParaMonte.spec";
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %function self = SpecMCMC()
        %end

        function result = chainSize(self,chainSize)
            result = verifySpec(chainSize,self.objectName+"chainSize","integer",delim);
        end

        function result = startPointVec(self,startPointVec)
            result = verifySpec(startPointVec,self.objectName+"startPointVec","real",delim);
        end

        function result = sampleRefinementCount(self,sampleRefinementCount)
            result = verifySpec(sampleRefinementCount,self.objectName+"sampleRefinementCount","integer",delim);
        end

        function result = sampleRefinementMethod(self,sampleRefinementMethod)
            result = verifySpec(sampleRefinementMethod,self.objectName+"sampleRefinementMethod","string",delim);
        end

        function result = randomStartPointRequested(self,randomStartPointRequested)
            result = verifySpec(randomStartPointRequested,self.objectName+"randomStartPointRequested","logical",delim);
        end

        function result = randomStartPointDomainLowerLimitVec(self,randomStartPointDomainLowerLimitVec)
            result = verifySpec(randomStartPointDomainLowerLimitVec,self.objectName+"randomStartPointDomainLowerLimitVec","real",delim);
        end

        function result = randomStartPointDomainUpperLimitVec(self,randomStartPointDomainUpperLimitVec)
            result = verifySpec(randomStartPointDomainUpperLimitVec,self.objectName+"randomStartPointDomainUpperLimitVec","real",delim);
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecMCMC