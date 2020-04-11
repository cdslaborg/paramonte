#**********************************************************************************************************************************
#**********************************************************************************************************************************
#
#  ParaMonte: plain powerful parallel Monte Carlo library.
#
#  Copyright (C) 2012-present, The Computational Data Science Lab
#
#  This file is part of ParaMonte library. 
#
#  ParaMonte is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, version 3 of the License.
#
#  ParaMonte is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
#
#**********************************************************************************************************************************
#**********************************************************************************************************************************

import os as _os
import numpy as _np
import sys as _sys
_sys.path.append(_os.path.dirname(__file__))

from _SpecBase import _delim

####################################################################################################################################
#### SpecDRAM specification type-checking class
####################################################################################################################################

class _SpecDRAM():

    def scaleFactor(self,scaleFactor):
        if isinstance(scaleFactor,str):
            return "scaleFactor=" + "'" + str(scaleFactor) + "'" + _delim
        else:
            raise TypeError("The input specification, scaleFactor, must be of type str.")

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

    def adaptiveUpdateCount(self,adaptiveUpdateCount):
        if isinstance(adaptiveUpdateCount,int):
            return "adaptiveUpdateCount=" + str(adaptiveUpdateCount) + _delim
        else:
            raise TypeError("The input specification, adaptiveUpdateCount, must be of type int.")

    def adaptiveUpdatePeriod(self,adaptiveUpdatePeriod):
        if isinstance(adaptiveUpdatePeriod,int):
            return "adaptiveUpdatePeriod=" + str(adaptiveUpdatePeriod) + _delim
        else:
            raise TypeError("The input specification, adaptiveUpdatePeriod, must be of type int.")

    def greedyAdaptationCount(self,greedyAdaptationCount):
        if isinstance(greedyAdaptationCount,int):
            return "greedyAdaptationCount=" + str(greedyAdaptationCount) + _delim
        else:
            raise TypeError("The input specification, greedyAdaptationCount, must be of type int.")

    def delayedRejectionCount(self,delayedRejectionCount):
        if isinstance(delayedRejectionCount,int):
            return "delayedRejectionCount=" + str(delayedRejectionCount) + _delim
        else:
            raise TypeError("The input specification, delayedRejectionCount, must be of type int.")

    def burninAdaptationMeasure(self,burninAdaptationMeasure):
        if isinstance(burninAdaptationMeasure,(int,float)):
            return "burninAdaptationMeasure=" + str(burninAdaptationMeasure) + _delim
        else:
            raise TypeError("The input specification, burninAdaptationMeasure, must be of type float.")

    def delayedRejectionScaleFactorVec(self,delayedRejectionScaleFactorVec):
        if isinstance(delayedRejectionScaleFactorVec,(list,tuple,_np.ndarray)):
            return "delayedRejectionScaleFactorVec=" + str(_np.array(list(delayedRejectionScaleFactorVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("The input specification, delayedRejectionScaleFactorVec must be list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################