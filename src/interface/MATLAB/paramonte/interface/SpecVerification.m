classdef SpecVerification

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    properties (Access = public)
        delim = ",";
        objectName = "";
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods(Access = public)

        function verifiedSpecString = verifySpec(self,variableValue,variableType)

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
                error("Internal Error: unrecognized variableType: " + string(variableType) );
            end

            objectNameCorrected = self.objectName;
            verifiedSpecString = variableName + "=";
            variableValueLen = length(variableValue);

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
                            error   ( "The input specification, " + objectNameCorrected + ", cannot contain both single-quote and double-quote characters. " ...
                                    + "Use only one type of quotation marks in your input string. " + objectNameCorrected + " = " ...
                                    + string(value)  ...
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
                    else
                        typeIsWrong = true;
                    end

                end

                if typeIsWrong
                    error("The input specification, " + objectNameCorrected + ", must be of type " + possibleTypes + ".");
                end

            end % for loop

        end % function

    end % methods

end % classdef SpecVerification
