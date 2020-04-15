classdef SpecBase_Description_class < handle

    properties (Constant)
        CLASS_NAME          = "@SpecBase_Description_class"
        MAX_DESCRIPTION_LEN = 4096
    end

    properties
        val                 = []
        def                 = []
        desc                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_Description_class(methodName)
            self.def    = "Nothing provided by the user.";
            self.desc   = "The variable 'description' contains general information about the specific " + methodName + " simulation that "              ...
                        + "is going to be performed. It has no effects on the simulation and serves only as a general description of the simulation "   ...
                        + "for future reference. The " + methodName + " parser automatically recognizes the C-style '\\n' escape sequence as the "      ...
                        + "new-line character, and '\\\\' as the backslash character '\\' if they used in the description. For example, '\\\\n' "       ...
                        + "will be converted to '\\n' on the output, while '\\n' translates to the new-line character. Other C escape sequences "       ...
                        + "are neither supported nor needed. The default value for description is '" + self.def + "'."                                  ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, description)
        
            if isempty(description)
                self.val = strtrim(self.def);
            else
                self.val = strtrim(description);
            end
            
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%**********************************************************************************************************************************

end