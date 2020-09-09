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
