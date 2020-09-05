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
#### SpecDRAM specification type-checking class
####################################################################################################################################

def adaptiveUpdateCount(self,adaptiveUpdateCount):
    if isinstance(adaptiveUpdateCount,int):
        return "adaptiveUpdateCount=" + str(adaptiveUpdateCount) + inputFileDelim
    else:
        raise TypeError("The input specification, adaptiveUpdateCount, must be of type int.")

####################################################################################################################################

def adaptiveUpdatePeriod(self,adaptiveUpdatePeriod):
    if isinstance(adaptiveUpdatePeriod,int):
        return "adaptiveUpdatePeriod=" + str(adaptiveUpdatePeriod) + inputFileDelim
    else:
        raise TypeError("The input specification, adaptiveUpdatePeriod, must be of type int.")

####################################################################################################################################

def greedyAdaptationCount(self,greedyAdaptationCount):
    if isinstance(greedyAdaptationCount,int):
        return "greedyAdaptationCount=" + str(greedyAdaptationCount) + inputFileDelim
    else:
        raise TypeError("The input specification, greedyAdaptationCount, must be of type int.")

####################################################################################################################################

def delayedRejectionCount(self,delayedRejectionCount):
    if isinstance(delayedRejectionCount,int):
        return "delayedRejectionCount=" + str(delayedRejectionCount) + inputFileDelim
    else:
        raise TypeError("The input specification, delayedRejectionCount, must be of type int.")

####################################################################################################################################

def burninAdaptationMeasure(self,burninAdaptationMeasure):
    if isinstance(burninAdaptationMeasure,(int,float)):
        return "burninAdaptationMeasure=" + str(burninAdaptationMeasure) + inputFileDelim
    else:
        raise TypeError("The input specification, burninAdaptationMeasure, must be of type float.")

####################################################################################################################################

def delayedRejectionScaleFactorVec(self,delayedRejectionScaleFactorVec):
    if isinstance(delayedRejectionScaleFactorVec,(list,tuple,np.ndarray)):
        return "delayedRejectionScaleFactorVec=" + str(np.array(list(delayedRejectionScaleFactorVec)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, delayedRejectionScaleFactorVec must be list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################

def helpme(specification = ""):
    """
    |
    Return help on the input ParaDRAM specification.
    If no input is provided, then the weblink to the 
    entire list of ParaDRAM specifications will be output.

        **Parameters**

            specification
                A string that can take be any of the simulation specifications 
                of a ParaDRAM object, such as, "chainSize", ...

        **Returns**

            None

    """
    import _paramonte as pm
    if isinstance(specification,str):
        if specification!="": specification = "#" + specification.lower()
    else:
        pm.abort( msg   = "The input argument to helpme() function must be a string whose value \n"
                        + "can be the name of one of the specifications of the ParaDRAM sampler.\n"
                        + "For example, to get information about the variable `chainSize`, try:\n\n"
                        + "    pmpd.spec.helpme(\"chainSize\")\n\n"
                        + "where ``pmpd`` is the name of the ParaDRAM sampler object.\n"
                        + "For more information and examples on the usage, visit:\n\n"
                        + "    " + pm.website.home.url
                , methodName = pm.names.paradram
                , marginTop = 1
                , marginBot = 1
                )

    print("Visit: " + pm.website.home.usage.paradram.specifications.url + specification)

####################################################################################################################################
