####################################################################################################################################
####################################################################################################################################
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   ParaMonte is free software: you can redistribute it and/or modify it
####   under the terms of the GNU Lesser General Public License as published
####   by the Free Software Foundation, version 3 of the License.
####
####   ParaMonte is distributed in the hope that it will be useful,
####   but WITHOUT ANY WARRANTY; without even the implied warranty of
####   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
####   GNU Lesser General Public License for more details.
####
####   You should have received a copy of the GNU Lesser General Public License
####   along with the ParaMonte library. If not, see,
####
####       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
####
####   ACKNOWLEDGMENT
####
####   As per the ParaMonte library license agreement terms,
####   if you use any parts of this library for any purposes,
####   we ask you to acknowledge the ParaMonte library's usage
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import os as _os
import numpy as _np
import sys as _sys
_sys.path.append(_os.path.dirname(__file__))

from _SpecBase import _delim

####################################################################################################################################
#### SpecMCMC specification type-checking class
####################################################################################################################################

class _SpecMCMC():

    def chainSize(self,chainSize):
        if isinstance(chainSize,int):
            return "chainSize=" + str(chainSize) + _delim
        else:
            raise TypeError("The input specification, chainSize, must be of type int.")

    def scaleFactor(self,scaleFactor):
        if isinstance(scaleFactor,str):
            return "scaleFactor=" + "'" + str(scaleFactor) + "'" + _delim
        else:
            raise TypeError("The input specification, scaleFactor, must be of type str.")

    def startPointVec(self,startPointVec):
        if isinstance(startPointVec,(list,tuple,_np.ndarray)):
            return "startPointVec=" + str(_np.array(list(startPointVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("The input specification, startPointVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

    def proposalModel(self,proposalModel):
        if isinstance(proposalModel,str):
            return "proposalModel=" + "'" + str(proposalModel) + "'" + _delim
        else:
            raise TypeError("The input specification, proposalModel, must be of type str.")

    def proposalStartCovMat(self,proposalStartCovMat):
        if isinstance(proposalStartCovMat,(list,tuple,_np.ndarray)):
            return "proposalStartCovMat=" + str(_np.array(list(proposalStartCovMat)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("The input specification, proposalStartCovMat, must be a symmetric 2-dimensional list, tuple, or numpy matrix of ndim-by-ndim or, "
                            "a list, tuple, or numpy vector of ndim**2 or less elements of type float.")

    def proposalStartCorMat(self,proposalStartCorMat):
        if isinstance(proposalStartCorMat,(list,tuple,_np.ndarray)):
            return "proposalStartCorMat=" + str(_np.array(list(proposalStartCorMat)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("The input specification, proposalStartCorMat, must be a symmetric 2-dimensional list, tuple, or numpy matrix of ndim-by-ndim or, "
                            "a list, tuple, or numpy vector of ndim**2 or less elements of type float.")

    def proposalStartStdVec(self,proposalStartStdVec):
        if isinstance(proposalStartStdVec,(list,tuple,_np.ndarray)):
            return "proposalStartStdVec=" + str(_np.array(list(proposalStartStdVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("The input specification, proposalStartStdVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

    def sampleRefinementCount(self,sampleRefinementCount):
        if isinstance(sampleRefinementCount,int):
            return "sampleRefinementCount=" + str(sampleRefinementCount) + _delim
        else:
            raise TypeError("The input specification, sampleRefinementCount, must be of type int.")

    def sampleRefinementMethod(self,sampleRefinementMethod):
        if isinstance(sampleRefinementMethod,str):
            return "sampleRefinementMethod=" + "'" + str(sampleRefinementMethod) + "'" + _delim
        else:
            raise TypeError("The input specification, sampleRefinementMethod, must be of type str.")

    def randomStartPointRequested(self,randomStartPointRequested):
        if isinstance(randomStartPointRequested,bool):
            return "randomStartPointRequested=" + str(randomStartPointRequested) + _delim
        else:
            raise TypeError("The input specification, randomStartPointRequested, must be of type bool (True or False).")

    def randomStartPointDomainLowerLimitVec(self,randomStartPointDomainLowerLimitVec):
        if isinstance(randomStartPointDomainLowerLimitVec,(list,tuple,_np.ndarray)):
            return "randomStartPointDomainLowerLimitVec=" + str(_np.array(list(randomStartPointDomainLowerLimitVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("The input specification, randomStartPointDomainLowerLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

    def randomStartPointDomainUpperLimitVec(self,randomStartPointDomainUpperLimitVec):
        if isinstance(randomStartPointDomainUpperLimitVec,(list,tuple,_np.ndarray)):
            return "randomStartPointDomainUpperLimitVec=" + str(_np.array(list(randomStartPointDomainUpperLimitVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("The input specification, randomStartPointDomainUpperLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################