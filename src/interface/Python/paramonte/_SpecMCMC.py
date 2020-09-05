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
####   we ask you to acknowledge the use of the ParaMonte library
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import numpy as np
from _SpecBase import inputFileDelim

####################################################################################################################################
#### SpecMCMC specification type-checking class
####################################################################################################################################

def chainSize(chainSize):
    if isinstance(chainSize,int):
        return "chainSize=" + str(chainSize) + inputFileDelim
    else:
        raise TypeError("The input specification, chainSize, must be of type int.")

####################################################################################################################################

def scaleFactor(scaleFactor):
    if isinstance(scaleFactor,str):
        return "scaleFactor=" + "'" + str(scaleFactor) + "'" + inputFileDelim
    else:
        raise TypeError("The input specification, scaleFactor, must be of type str.")

####################################################################################################################################

def startPointVec(startPointVec):
    if isinstance(startPointVec,(list,tuple,np.ndarray)):
        return "startPointVec=" + str(np.array(list(startPointVec)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, startPointVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################

def proposalModel(proposalModel):
    if isinstance(proposalModel,str):
        return "proposalModel=" + "'" + str(proposalModel) + "'" + inputFileDelim
    else:
        raise TypeError("The input specification, proposalModel, must be of type str.")

####################################################################################################################################

def proposalStartCovMat(proposalStartCovMat):
    if isinstance(proposalStartCovMat,(list,tuple,np.ndarray)):
        return "proposalStartCovMat=" + str(np.array(list(proposalStartCovMat)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, proposalStartCovMat, must be a symmetric 2-dimensional list, tuple, or numpy matrix of ndim-by-ndim or, "
                        "a list, tuple, or numpy vector of ndim**2 or less elements of type float.")

####################################################################################################################################

def proposalStartCorMat(proposalStartCorMat):
    if isinstance(proposalStartCorMat,(list,tuple,np.ndarray)):
        return "proposalStartCorMat=" + str(np.array(list(proposalStartCorMat)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, proposalStartCorMat, must be a symmetric 2-dimensional list, tuple, or numpy matrix of ndim-by-ndim or, "
                        "a list, tuple, or numpy vector of ndim**2 or less elements of type float.")

####################################################################################################################################

def proposalStartStdVec(proposalStartStdVec):
    if isinstance(proposalStartStdVec,(list,tuple,np.ndarray)):
        return "proposalStartStdVec=" + str(np.array(list(proposalStartStdVec)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, proposalStartStdVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################

def sampleRefinementCount(sampleRefinementCount):
    if isinstance(sampleRefinementCount,int):
        return "sampleRefinementCount=" + str(sampleRefinementCount) + inputFileDelim
    else:
        raise TypeError("The input specification, sampleRefinementCount, must be of type int.")

####################################################################################################################################

def sampleRefinementMethod(sampleRefinementMethod):
    if isinstance(sampleRefinementMethod,str):
        return "sampleRefinementMethod=" + "'" + str(sampleRefinementMethod) + "'" + inputFileDelim
    else:
        raise TypeError("The input specification, sampleRefinementMethod, must be of type str.")

####################################################################################################################################

def randomStartPointRequested(randomStartPointRequested):
    if isinstance(randomStartPointRequested,bool):
        return "randomStartPointRequested=" + str(randomStartPointRequested) + inputFileDelim
    else:
        raise TypeError("The input specification, randomStartPointRequested, must be of type bool (True or False).")

####################################################################################################################################

def randomStartPointDomainLowerLimitVec(randomStartPointDomainLowerLimitVec):
    if isinstance(randomStartPointDomainLowerLimitVec,(list,tuple,np.ndarray)):
        return "randomStartPointDomainLowerLimitVec=" + str(np.array(list(randomStartPointDomainLowerLimitVec)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, randomStartPointDomainLowerLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################

def randomStartPointDomainUpperLimitVec(randomStartPointDomainUpperLimitVec):
    if isinstance(randomStartPointDomainUpperLimitVec,(list,tuple,np.ndarray)):
        return "randomStartPointDomainUpperLimitVec=" + str(np.array(list(randomStartPointDomainUpperLimitVec)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, randomStartPointDomainUpperLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################

