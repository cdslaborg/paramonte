classdef Err_class < handle

    properties (Constant)
        CLASS_NAME              = "@Err_class"
        ERR_HANDLING_REQUESTED  = false
    end

    properties
        occurred                = []
        stat                    = []
        statNull                = []
        msg                     = ''
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = Err_class()
            self.occurred    = false;
            self.stat        = -intmax;
            self.statNull    = -intmax;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function abort(self, prefix, newLine, outputUnit)

            Decoration  = Decoration_class([],[],[],[]);
            imageChar   = num2str(1);
            msg         = self.msg;

            if ~isempty(prefix)
                self.informUser(msg, prefix + " - FATAL: "  , newLine, outputUnit, [], [], [], []);
                pfx = prefix;
            else
                self.informUser(msg, " - "                  , newLine, outputUnit, [], [], [], []);
                pfx = "";
            end

            if ~isempty(outputUnit)
                if outputUnit ~= 1
                    Decoration.write(outputUnit, 1, 0, 1, pfx + " - Please Correct the error(s) and rerun the simulation.");
                    Decoration.write(outputUnit, 1, 0, 1, pfx + " - For further help, contact Amir Shahmoradi via:");
                    Decoration.write(outputUnit, 0, 0, 1, pfx + " - a.shahmoradi@gmail.com");
                    Decoration.write(outputUnit, 0, 0, 1, pfx + " - shahmoradi@utexas.edu");
                    Decoration.write(outputUnit, 0, 0, 1, pfx + " - cdslab.org/ParaMonte/");
                    Decoration.write(outputUnit, 1, 2, 1, pfx + " - Gracefully Exiting on image " + strtrim(imageChar) + ".");
                end
            end

            if outputUnit ~= 1
                % notify the user on screen too
                Decoration.write(1, 1, 0, 1, pfx + " - FATAL: Runtime error occurred.");
                Decoration.write(1, 0, 0, 1, pfx + " - FATAL: For more information please see the report file.");
                Decoration.write(1, 0, 2, 1, pfx + " - FATAL: Gracefully Exiting on image " + strtrim(imageChar) + ".");
            end
            
            error stop
            
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods(Static)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function warn(msg, prefix, newLine, outputUnit)
            if ~isempty(prefix)
                Err_class.informUser(msg, prefix + " - WARNING: " , newLine, outputUnit, [], [], [], []);
            else
                Err_class.informUser(msg, " - WARNING: "          , newLine, outputUnit, [], [], [], []);
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function note(msg, prefix, newLine, outputUnit, marginTop, marginBot)

            if ~isempty(prefix)
                Err_class.informUser(msg, prefix + " - NOTE: ", newLine, outputUnit, [], [], marginTop, marginBot);
            else
                Err_class.informUser(msg, " - NOTE: "         , newLine, outputUnit, [], [], marginTop, marginBot);
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function informUser(msg, prefix, newLine, outputUnit, wrapSplit, wrapWidth, marginTop, marginBot)
            
            msg         = strrep(msg, '\n', newline);
            Decoration  = Decoration_class([],[],[],[]);

            if ~isempty(outputUnit) ,   stdout  = outputUnit;   else,   stdout  = 1     ; end
            if ~isempty(prefix)     ,   pfx     = prefix    ;   else,   pfx     = ""    ; end
            if ~isempty(wrapSplit)  ,   split   = wrapSplit ;   else,   split   = " "   ; end
            if ~isempty(wrapWidth)  ,   width   = wrapWidth ;   else,   width   = 100   ; end
            if ~isempty(marginTop)  ,   padTop  = marginTop ;   else,   padTop  = 1     ; end
            if ~isempty(marginBot)  ,   padBot  = marginBot ;   else,   padBot  = 1     ; end

            List = Decoration.getListOfLines(msg, newline);

            lenList = length(List);
            for i = 1 : lenList
                ListJustified       = Decoration.wrapText(List{i}, width, split , []);
                lenListJustified    = length(ListJustified);
                for ijustified = 1 : lenListJustified
                    padTopCurrent = 0;
                    padBotCurrent = 0;
                    if (i == 1) && (ijustified == 1)                        , padTopCurrent = padTop; end % the very first line
                    if (i == lenList) && (ijustified == lenListJustified)   , padBotCurrent = padBot; end % the very last line
                    Decoration.write(stdout, padTopCurrent, padBotCurrent   , 1, pfx + ListJustified{ijustified});
                end
            end

            if isempty(marginBot), Decoration.write(stdout,[],[],[],[]); end

        end % function informUser

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end