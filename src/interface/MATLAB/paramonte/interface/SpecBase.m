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
%%%% SpecBase specification type-checking class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef SpecBase

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        delim = ",";
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %function self = SpecBase()
        %end

        function result = sampleSize(self,sampleSize)
            if isa(sampleSize,"numeric")
                result = "sampleSize=" + num2str(sampleSize) + self.delim
            else
                error("The input specification, sampleSize, must be numeric.")
            end
        end

        function result = randomSeed(self,randomSeed)
            if isa(randomSeed,"numeric")
                result = "randomSeed=" + num2str(int32(randomSeed)) + self.delim
            else
                error("The input specification, randomSeed, must be numeric.")
            end
        end

            else:
                raise TypeError("The input specification, description, must be of type str.")

        function result = description(self,description)
            if isa(description,'string') || isa(description,'char')
                description = string(description)
                hasSingleQuote = contains(description,"'");
                hasDoubleQuote = contains(description,'"');
                if hasSingleQuote
                    if hasDoubleQuote
                        error   ( "The input specification, description, cannot contain both single-quote and double-quote characters. " ...
                                + "Use only one type of quotation marks in your input string. description = " ...
                                + description  ...
                                )
                    else
                        result = "description=" + """" + description + """" + _delim
                    end
                else
                    result = "description=" + "'" + description + "'" + _delim
                end
                result = "description=" + "'" + description + "'" + self.delim
            else
                error("The input specification, description, must be of type string or char.")
            end
        end

        function result = outputFileName(self,outputFileName)
            if isa(outputFileName,'string') || isa(outputFileName,'char')
                result = "outputFileName=" + "'" + outputFileName + "'" + self.delim
            else
                error("The input specification, outputFileName, must be of type str.")
            end
        end

        function result = outputDelimiter(self,outputDelimiter)
            if isa(outputDelimiter,str)
                result = "outputDelimiter=" + "'" + num2str(outputDelimiter) + "'" + self.delim
            else
                error("The input specification, outputDelimiter, must be of type str.")
            end
        end

        function result = chainFileFormat(self,chainFileFormat)
            if isa(chainFileFormat,str)
                result = "chainFileFormat=" + "'" + num2str(chainFileFormat) + "'" + self.delim
            else
                error("The input specification, chainFileFormat, must be of type str.")
            end
        end

        function result = variableNameList(self,variableNameList)
            if isa(variableNameList,(list,tuple))
                result = "variableNameList=" + self.delim.join("'{0}'".format(_) for _ in list(variableNameList)) + self.delim
            else
                error("The input specification, variableNameList, must be either a list or tuple of ndim or less elements, each element of which must be of type str.")
            end
        end

        function result = restartFileFormat(self,restartFileFormat)
            if isa(restartFileFormat,str)
                result = "restartFileFormat=" + "'" + num2str(restartFileFormat) + "'" + self.delim
            else
                error("The input specification, restartFileFormat, must be of type str.")
            end
        end

        function result = outputColumnWidth(self,outputColumnWidth)
            if isa(outputColumnWidth,int)
                result = "outputColumnWidth=" + num2str(outputColumnWidth) + self.delim
            else
                error("The input specification, outputColumnWidth, must be of type int.")
            end
        end

        function result = outputRealPrecision(self,outputRealPrecision)
            if isa(outputRealPrecision,int)
                result = "outputRealPrecision=" + num2str(outputRealPrecision) + self.delim
            else
                error("The input specification, outputRealPrecision, must be of type int.")
            end
        end

        function result = silentModeRequested(self,silentModeRequested)
            if isa(silentModeRequested,bool)
                result = "silentModeRequested=" + num2str(silentModeRequested) + self.delim
            else
                error("The input specification, silentModeRequested, must be of type bool (True or False).")
            end
        end

        function result = domainLowerLimitVec(self,domainLowerLimitVec)
            if isa(domainLowerLimitVec,(list,tuple,_np.ndarray))
                result = "domainLowerLimitVec=" + num2str(_np.array(list(domainLowerLimitVec)).flatten()).strip('[]') + self.delim
            else
                error("The input specification, domainLowerLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")
            end
        end

        function result = domainUpperLimitVec(self,domainUpperLimitVec)
            if isa(domainUpperLimitVec,(list,tuple,_np.ndarray))
                result = "domainUpperLimitVec=" + num2str(_np.array(list(domainUpperLimitVec)).flatten()).strip('[]') + self.delim
            else
                error("The input specification, domainUpperLimitVec, must be a list, tuple, or numpy vector of ndim or less elements of type float.")
            end
        end

        function result = parallelizationModel(self,parallelizationModel)
            if isa(parallelizationModel,str)
                result = "parallelizationModel=" + "'" + num2str(parallelizationModel) + "'" + self.delim
            else
                error("The input specification, parallelizationModel, must be of type str.")
            end
        end

        function result = progressReportPeriod(self,progressReportPeriod)
            if isa(progressReportPeriod,int)
                result = "progressReportPeriod=" + num2str(progressReportPeriod) + self.delim
            else
                error("The input specification, progressReportPeriod, must be of type int.")
            end
        end

        function result = targetAcceptanceRate(self,targetAcceptanceRate)
            if isa(targetAcceptanceRate,float)
                result = "targetAcceptanceRate=" + num2str(targetAcceptanceRate) + self.delim
            else
                error("The input specification, targetAcceptanceRate, must be of type float.")
            end
        end

        function result = mpiFinalizeRequested(self,mpiFinalizeRequested)
            if isa(mpiFinalizeRequested,bool)
                result = "mpiFinalizeRequested=" + num2str(mpiFinalizeRequested) + self.delim
            else
                error("The input specification, mpiFinalizeRequested, must be of type bool (True or False).")
            end
        end

        function result = maxNumDomainCheckToWarn(self,maxNumDomainCheckToWarn)
            if isa(maxNumDomainCheckToWarn,int)
                result = "maxNumDomainCheckToWarn=" + num2str(maxNumDomainCheckToWarn) + self.delim
            else
                error("The input specification, maxNumDomainCheckToWarn, must be of type int.")
            end
        end

        function result = maxNumDomainCheckToStop(self,maxNumDomainCheckToStop)
            if isa(maxNumDomainCheckToStop,int)
                result = "maxNumDomainCheckToStop=" + num2str(maxNumDomainCheckToStop) + self.delim
            else
                error("maxNumDomainCheckToStop must be of type int.")
            end
        end

        function result = interfaceType(self)
            result = "interfaceType=" + "'Python " + _sys.version + "'" + self.delim
        end

    end % methods (dynamic)

####################################################################################################################################
#### generate output filename
####################################################################################################################################

function result = _genOutputFileName(methodName)
    from datetime import datetime as _dt
    dt = _dt.now()
    result  methodName + "_run_" \
            + "{:04d}".format(dt.year) + "{:02d}".format(dt.month) + "{:02d}".format(dt.day) + "_"  \
            + "{:02d}".format(dt.hour) + "{:02d}".format(dt.minute) + "{:02d}".format(dt.second) + "_" \
            + "{:03d}".format(round(dt.microsecond/1000))

####################################################################################################################################
