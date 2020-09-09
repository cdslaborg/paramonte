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