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
#### SpecMCMC specification type-checking class
####################################################################################################################################

class _SpecMCMC():

    def chainSize(self,chainSize):
        if isinstance(chainSize,int):
            return "chainSize=" + str(chainSize) + _delim
        else:
            raise TypeError("chainSize must be of type int.")

    def startPointVec(self,startPointVec):
        if isinstance(startPointVec,(list,tuple,_np.ndarray)):
            return "startPointVec=" + str(_np.array(list(startPointVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("startPointVec must be a list, tuple, or numpy vector of ndim or less elements of type float.")

    def sampleRefinementCount(self,sampleRefinementCount):
        if isinstance(sampleRefinementCount,int):
            return "sampleRefinementCount=" + str(sampleRefinementCount) + _delim
        else:
            raise TypeError("sampleRefinementCount must be of type int.")

    def sampleRefinementMethod(self,sampleRefinementMethod):
        if isinstance(sampleRefinementMethod,str):
            return "sampleRefinementMethod=" + "'" + str(sampleRefinementMethod) + "'" + _delim
        else:
            raise TypeError("sampleRefinementMethod must be of type str.")

    def randomStartPointRequested(self,randomStartPointRequested):
        if isinstance(randomStartPointRequested,bool):
            return "randomStartPointRequested=" + str(randomStartPointRequested) + _delim
        else:
            raise TypeError("randomStartPointRequested must be of type bool (True or False).")

    def randomStartPointDomainLowerLimitVec(self,randomStartPointDomainLowerLimitVec):
        if isinstance(randomStartPointDomainLowerLimitVec,(list,tuple,_np.ndarray)):
            return "randomStartPointDomainLowerLimitVec=" + str(_np.array(list(randomStartPointDomainLowerLimitVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("randomStartPointDomainLowerLimitVec must be a list, tuple, or numpy vector of ndim or less elements of type float.")

    def randomStartPointDomainUpperLimitVec(self,randomStartPointDomainUpperLimitVec):
        if isinstance(randomStartPointDomainUpperLimitVec,(list,tuple,_np.ndarray)):
            return "randomStartPointDomainUpperLimitVec=" + str(_np.array(list(randomStartPointDomainUpperLimitVec)).flatten()).strip('[]') + _delim
        else:
            raise TypeError("randomStartPointDomainUpperLimitVec must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################