classdef SpecBase_ProgressReportPeriod_class < handle

    properties (Constant)
        CLASS_NAME  = "@SpecBase_ProgressReportPeriod_class"
    end

    properties
        val         = []
        def         = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_ProgressReportPeriod_class()
            self.def 	= 1000;
            self.desc	= "Every progressReportPeriod calls to the objective function, the sampling progress will be reported to the log file. "...
						+ "Note that progressReportPeriod must be a positive integer. "                                                         ...
						+ "The default value is " + num2str(self.def) + "."                                                                     ...
						;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, progressReportPeriod)
            self.val = progressReportPeriod;
            if isempty(self.val)
                self.val = self.def;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";
            if self.val < 1
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                                        ...
                                + "The input value for variable progressReportPeriod must be a positive integer value. If you are not sure "    ...
                                + "about the appropriate value for this variable, simply drop it from the input. "                              ...
                                + methodName + " will automatically assign an appropriate value to it." + newline + newline                     ...
                                ;
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end