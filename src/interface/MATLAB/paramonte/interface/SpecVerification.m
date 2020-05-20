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
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef SpecVerification

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        ndim = [];
        delim = ",";
        objectName = "";
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods(Access = public)

        function verifiedSpecString = verifySpec(self,variableValue,variableType,maxVariableLength)

            % verify string value of the input variable does not contain both single and double quotations.

            variableName = inputname(2);

            typeMustBeReal = false;
            typeMustBeString = false;
            typeMustBeInteger = false;
            typeMustBeLogical = false;

            if strcmp(variableType,"real")
                typeMustBeReal = true;
                possibleTypes = "numeric";
            elseif strcmp(variableType,"string")
                typeMustBeString = true;
                possibleTypes = "string or char";
                variableValue = string(variableValue);
            elseif strcmp(variableType,"integer")
                typeMustBeInteger = true;
                possibleTypes = "numeric";
            elseif strcmp(variableType,"logical")
                typeMustBeLogical = true;
                possibleTypes = "logical";
            else
                error(newline + "Internal Error: unrecognized variableType: " + string(variableType) + newline + newline );
            end

            objectNameCorrected = self.objectName + ".spec." + variableName;
            verifiedSpecString = variableName + "=";
            variableValueLen = length(variableValue);

            if variableValueLen > maxVariableLength

                variableValueList = "";
                for i = 1:variableValueLen
                    variableValueList = variableValueList + string(variableValue(i));
                    if i<variableValueLen; variableValueList = variableValueList + " "; end
                end
                error   ( newline ...
                        + "The length of the value of input specification, " + objectNameCorrected + ", cannot be larger than " + string(maxVariableLength) + "." + newline ...
                        + "You have specified:" + newline + newline ...
                        + "    " + objectNameCorrected + " = """ + string(variableValueList) + """" + newline ...
                        );

            else

                for i = 1:variableValueLen

                    if variableValueLen>1; objectNameCorrected = self.objectName + "(" + string(i) + ")"; end;

                    if isa(variableValue(i),"cell")
                        value = variableValue{i};
                    else
                        value = variableValue(i);
                    end

                    typeIsWrong = false;

                    if typeMustBeReal

                        if isa(value,"numeric")
                            verifiedSpecString = verifiedSpecString + value + self.delim;
                        else
                            typeIsWrong = true;
                        end

                    elseif typeMustBeString

                        if isa(value,"string") || isa(value,"char")
                            enclosedString = encloseString( string(value) );
                            if isempty(enclosedString)
                                error   ( newline ...
                                        + "The input specification, " + objectNameCorrected + ", cannot contain both single-quote and double-quote characters. " ...
                                        + "Use only one type of quotation marks in your input string. " + objectNameCorrected + " = " ...
                                        + string(value)  ...
                                        + newline+ newline ...
                                        );
                            else
                                verifiedSpecString = verifiedSpecString + enclosedString + self.delim;
                            end
                        else
                            typeIsWrong = true;
                        end

                    elseif typeMustBeInteger

                        if isa(value,"numeric")
                            verifiedSpecString = verifiedSpecString + int32(value) + self.delim;
                        else
                            typeIsWrong = true;
                        end

                    elseif typeMustBeLogical

                        if isa(value,"logical")
                            verifiedSpecString = verifiedSpecString + string(value) + self.delim;
                        elseif isa(value,"numeric")
                            if value==0
                                verifiedSpecString = verifiedSpecString + ".false." + self.delim;
                            elseif value==1
                                verifiedSpecString = verifiedSpecString + ".true." + self.delim;
                            else
                                typeIsWrong = true;
                            end
                        else
                            typeIsWrong = true;
                        end

                    end

                    if typeIsWrong
                        error   ( newline ...
                                + "The input specification, " + objectNameCorrected + ", must be of type " + possibleTypes + ". You have entered:" + newline  + newline ...
                                + "    " + objectNameCorrected + " = " + string(value) + newline ...
                                + "    class(" + objectNameCorrected + ") = " + string(class(value)) + newline ...
                                );
                    end

                end % for loop

            end % if-block

        end % function

    end % methods

end % classdef SpecVerification
