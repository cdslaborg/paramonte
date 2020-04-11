classdef SpecBase_OutputFileName_class < Path_class

    properties (Constant)
        SUB_CLASS_NAME = "@SpecBase_OutputFileName_class"
    end

    properties
        def         = []
        namePrefix  = []
        pathPrefix  = []
        desc        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_OutputFileName_class(methodName)
            
            self = self@Path_class([], []);
            
            self.def    = methodName + '_run_' + datestr(now, 'ddmmyy') + '_' + datestr(now, 'HHMMSS_FFF');
            self.desc   = "outputFileName contains the path and the base of the filename for " + methodName + " output files. "                         ...
                        + "If not provided by the user, the default outputFileName is constructed from the current date and time:"                      ...
                        + newline + newline + Decoration_class.TAB + methodName + "_run_ddmmyy_hhmmss_fff" + newline + newline                          ...
                        + "where dd, mm, yy, hh, mm, ss stand respectively for the current day, month, year, hour, minute, second. "                    ...
                        + "In such a case, the default directory for the output files will be the current working directory of "                        ...
                        + methodName + ". If outputFileName is provided, but ends with a separator character '/' or '\' (as in Linux or Windows OS), "  ...
                        + "then its value will be used as the directory to which " + methodName + " output files will be written. In this case, "       ...
                        + "the output file naming convention described above will be used. Also, the given directory will be automatically created "    ...
                        + "if it does not exist already."                                                                                               ...
                        ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, outputFileName)
            if isempty(outputFileName)
                self.original = self.def;
            else
                self.original = strtrim(outputFileName);
            end
            self.modified = self.original;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end