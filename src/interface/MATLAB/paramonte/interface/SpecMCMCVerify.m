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
%   we ask you to acknowledge the ParaMonte library's usage
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SpecMCMC specification type-checking class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SpecMCMCVerify < SpecVerification

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        function self = SpecMCMCVerify(objectName,ndim)
            self.objectName = objectName;
            self.ndim = ndim;
        end

        function result = chainSize(self,chainSize)
            result = self.verifySpec(chainSize,"integer",1);
        end

        function result = scaleFactor(self,scaleFactor)
            result = self.verifySpec(scaleFactor,"string",1);
        end

        function result = startPointVec(self,startPointVec)
            result = self.verifySpec(startPointVec,"real",self.ndim);
        end

        function result = proposalModel(self,proposalModel)
            result = self.verifySpec(proposalModel,"string",1);
        end

        function result = proposalStartCovMat(self,proposalStartCovMat)
            result = self.verifySpec(proposalStartCovMat,"real",self.ndim^2);
        end

        function result = proposalStartCorMat(self,proposalStartCorMat)
            result = self.verifySpec(proposalStartCorMat,"real",self.ndim^2);
        end

        function result = proposalStartStdVec(self,proposalStartStdVec)
            result = self.verifySpec(proposalStartStdVec,"real",self.ndim);
        end

        function result = sampleRefinementCount(self,sampleRefinementCount)
            result = self.verifySpec(sampleRefinementCount,"integer",1);
        end

        function result = sampleRefinementMethod(self,sampleRefinementMethod)
            result = self.verifySpec(sampleRefinementMethod,"string",1);
        end

        function result = randomStartPointRequested(self,randomStartPointRequested)
            result = self.verifySpec(randomStartPointRequested,"logical",1);
        end

        function result = randomStartPointDomainLowerLimitVec(self,randomStartPointDomainLowerLimitVec)
            result = self.verifySpec(randomStartPointDomainLowerLimitVec,"real",self.ndim);
        end

        function result = randomStartPointDomainUpperLimitVec(self,randomStartPointDomainUpperLimitVec)
            result = self.verifySpec(randomStartPointDomainUpperLimitVec,"real",self.ndim);
        end

    end % methods (dynamic)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end % classdef SpecDRAMVerify