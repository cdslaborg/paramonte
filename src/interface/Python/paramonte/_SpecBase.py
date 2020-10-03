####################################################################################################################################
####################################################################################################################################
####
####   MIT License
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   Permission is hereby granted, free of charge, to any person obtaining a 
####   copy of this software and associated documentation files (the "Software"), 
####   to deal in the Software without restriction, including without limitation 
####   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
####   and/or sell copies of the Software, and to permit persons to whom the 
####   Software is furnished to do so, subject to the following conditions:
####
####   The above copyright notice and this permission notice shall be 
####   included in all copies or substantial portions of the Software.
####
####   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
####   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
####   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
####   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
####   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
####   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
####   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
####
####   ACKNOWLEDGMENT
####
####   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
####   As per the ParaMonte library license agreement terms, if you use any parts of 
####   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
####   work (education/research/industry/development/...) by citing the ParaMonte 
####   library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import sys
import numpy as np
from _pmutils import getList, getRandomFilePrefix

####################################################################################################################################
#### SpecBase specification type-checking class
####################################################################################################################################

inputFileDelim = ","

####################################################################################################################################

def sampleSize(sampleSize):
    if isinstance(sampleSize,(float,int)):
        return "sampleSize=" + str(sampleSize) + inputFileDelim
    else:
        raise TypeError("The input specification, sampleSize, must be of type float or int.")

####################################################################################################################################

def randomSeed(randomSeed):
    if isinstance(randomSeed,int):
        return "randomSeed=" + str(randomSeed) + inputFileDelim
    else:
        raise TypeError("The input specification, randomSeed, must be of type int.")

####################################################################################################################################

def description(description): return "description=" + verifyEncloseString(description,"description") + inputFileDelim

####################################################################################################################################

def outputFileName(outputFileName): return "outputFileName=" + verifyEncloseString(outputFileName,"outputFileName") + inputFileDelim

####################################################################################################################################

def outputDelimiter(outputDelimiter): return "outputDelimiter=" + verifyEncloseString(outputDelimiter,"outputDelimiter") + inputFileDelim

####################################################################################################################################

def chainFileFormat(chainFileFormat):  return "chainFileFormat=" + verifyEncloseString(chainFileFormat,"chainFileFormat") + inputFileDelim

####################################################################################################################################

def variableNameList(variableNameList):
    if isinstance(variableNameList,(list,tuple)):
        return "variableNameList=" + inputFileDelim.join("{0}".format(verifyEncloseString(_,"variableNameList["+str(i)+"]")) for i,_ in enumerate(getList(variableNameList))) + inputFileDelim
    else:
        raise TypeError("The input specification, variableNameList, must be either a list or tuple of ndim or less elements, each element of which must be of type str.")

####################################################################################################################################

def restartFileFormat(restartFileFormat): return "restartFileFormat=" + verifyEncloseString(restartFileFormat,"restartFileFormat") + inputFileDelim

####################################################################################################################################

def outputColumnWidth(outputColumnWidth):
    if isinstance(outputColumnWidth,int):
        return "outputColumnWidth=" + str(outputColumnWidth) + inputFileDelim
    else:
        raise TypeError("The input specification, outputColumnWidth, must be of type int.")

####################################################################################################################################

def overwriteRequested(overwriteRequested):
    if isinstance(overwriteRequested,bool):
        return "overwriteRequested=" + str(overwriteRequested) + inputFileDelim
    else:
        raise TypeError("The input specification, overwriteRequested, must be of type bool (True or False).")

####################################################################################################################################

def outputRealPrecision(outputRealPrecision):
    if isinstance(outputRealPrecision,int):
        return "outputRealPrecision=" + str(outputRealPrecision) + inputFileDelim
    else:
        raise TypeError("The input specification, outputRealPrecision, must be of type int.")

####################################################################################################################################

def silentModeRequested(silentModeRequested):
    if isinstance(silentModeRequested,bool):
        return "silentModeRequested=" + str(silentModeRequested) + inputFileDelim
    else:
        raise TypeError("The input specification, silentModeRequested, must be of type bool (True or False).")

####################################################################################################################################

def domainLowerLimitVec(domainLowerLimitVec):
    if isinstance(domainLowerLimitVec,(list,tuple,np.ndarray)):
        return "domainLowerLimitVec=" + str(np.array(getList(domainLowerLimitVec)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, domainLowerLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################

def domainUpperLimitVec(domainUpperLimitVec):
    if isinstance(domainUpperLimitVec,(list,tuple,np.ndarray)):
        return "domainUpperLimitVec=" + str(np.array(getList(domainUpperLimitVec)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, domainUpperLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")

####################################################################################################################################

def parallelizationModel(parallelizationModel):
    if isinstance(parallelizationModel,str):
        return "parallelizationModel=" + "'" + str(parallelizationModel) + "'" + inputFileDelim
    else:
        raise TypeError("The input specification, parallelizationModel, must be of type str.")

####################################################################################################################################

def progressReportPeriod(progressReportPeriod):
    if isinstance(progressReportPeriod,int):
        return "progressReportPeriod=" + str(progressReportPeriod) + inputFileDelim
    else:
        raise TypeError("The input specification, progressReportPeriod, must be of type int.")

####################################################################################################################################

def targetAcceptanceRate(targetAcceptanceRate):
    if isinstance(targetAcceptanceRate,(float,list,tuple,np.ndarray)):
        return "targetAcceptanceRate=" + str(np.array(getList(targetAcceptanceRate)).flatten()).strip('[]') + inputFileDelim
    else:
        raise TypeError("The input specification, targetAcceptanceRate, must be of type float.")

####################################################################################################################################

def mpiFinalizeRequested(mpiFinalizeRequested):
    if isinstance(mpiFinalizeRequested,bool):
        return "mpiFinalizeRequested=" + str(mpiFinalizeRequested) + inputFileDelim
    else:
        raise TypeError("The input specification, mpiFinalizeRequested, must be of type bool (True or False).")

####################################################################################################################################

def maxNumDomainCheckToWarn(maxNumDomainCheckToWarn):
    if isinstance(maxNumDomainCheckToWarn,int):
        return "maxNumDomainCheckToWarn=" + str(maxNumDomainCheckToWarn) + inputFileDelim
    else:
        raise TypeError("The input specification, maxNumDomainCheckToWarn, must be of type int.")

####################################################################################################################################

def maxNumDomainCheckToStop(maxNumDomainCheckToStop):
    if isinstance(maxNumDomainCheckToStop,int):
        return "maxNumDomainCheckToStop=" + str(maxNumDomainCheckToStop) + inputFileDelim
    else:
        raise TypeError("The input specification, maxNumDomainCheckToStop, must be of type int.")

####################################################################################################################################

def interfaceType():
    return "interfaceType=" + "'Python " + sys.version + "'" + inputFileDelim

####################################################################################################################################

def systemInfoFilePath(systemInfoFilePath):
    if systemInfoFilePath is None:
        return ""
    else:
        return "systemInfoFilePath='" + systemInfoFilePath + "'" + inputFileDelim

####################################################################################################################################
#### generate output filename
####################################################################################################################################

def genOutputFileName(methodName): return getRandomFilePrefix(prefix = methodName + "_run_")

####################################################################################################################################
#### enclose string with either single or double quotation . Return None if string contains both quotation symbols
####################################################################################################################################

def encloseString(string):
    hasSingleQuote = "'" in string
    hasDoubleQuote = '"' in string
    if hasSingleQuote:
        if hasDoubleQuote:
            return None
        else:
            return '"' + string + '"'
    else:
        return "'" + string + "'"

####################################################################################################################################
#### verify string value of the input variable does not contain both single and double quotations.
####################################################################################################################################

def verifyEncloseString(variableValue,variableName):
    if isinstance(variableValue,str):
        enclosedString = encloseString(variableValue)
        if enclosedString:
            return enclosedString
        else:
            raise TypeError ( "The input specification, " + variableName + ", cannot contain both single-quote and double-quote characters. "
                            + "Use only one type of quotation marks in your input string. " + variableName + " = "
                            + variableValue 
                            )
    else:
        raise TypeError("The input specification, " + variableName + ", must be of type str.")
