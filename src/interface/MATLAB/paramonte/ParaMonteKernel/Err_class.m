classdef Err_class < handle

    properties (Constant)
        CLASS_NAME              = "@Err_class"
        ERR_HANDLING_REQUESTED  = false
    end

    properties
        msg
        prefix
        wrapSplit
        wrapWidth
        marginTop
        marginBot
        outputUnit
        occurred
        newLine
        stat
        statNull
        resetEnabled
    end

    properties(Hidden)
        fullprefix
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Access = public)

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function self = Err_class()
            self.resetEnabled = true;
            self.reset();
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function abort(self)

            Decoration  = Decoration_class([],[],[],[]);
            imageChar   = num2str(1);

            if isempty(self.prefix)
                self.fullprefix = " - FATAL: ";
                self.informUser();
            else
                self.fullprefix = self.prefix + " - FATAL: ";
                self.informUser();
            end

            if self.outputUnit == 1
                % notify the user on screen too
                Decoration.write(1, 1, 0, 1, self.fullprefix + "Runtime error occurred.");
                Decoration.write(1, 0, 0, 1, self.fullprefix + "For more information please see the report file.");
                Decoration.write(1, 0, 2, 1, self.fullprefix + "Gracefully Exiting on image " + strtrim(imageChar) + ".");
            else
                Decoration.write(self.outputUnit, 1, 0, 1, self.fullprefix + "Please Correct the error(s) and rerun the simulation.");
                Decoration.write(self.outputUnit, 1, 0, 1, self.fullprefix + "For further help, contact the ParaMonte authors via:");
                Decoration.write(self.outputUnit, 0, 0, 1, self.fullprefix + "cdslab.org/ParaMonte/");
                Decoration.write(self.outputUnit, 1, 2, 1, self.fullprefix + "Gracefully Exiting on image " + strtrim(imageChar) + ".");
            end

            error stop

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function warn(self)
            if isempty(convertStringsToChars(self.prefix))
                self.fullprefix = " - WARNING: ";
                self.informUser();
            else
                self.fullprefix = self.prefix + " - WARNING: ";
                self.informUser();
            end
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function note(self)
            if isempty(convertStringsToChars(self.prefix))
                self.fullprefix = " - NOTE: ";
                self.informUser();
            else
                self.fullprefix = self.prefix + " - NOTE: ";
                self.informUser();
            end
        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        %function informUser(msg, prefix, newLine, outputUnit, wrapSplit, wrapWidth, marginTop, marginBot)
        function informUser(self)

            Decoration  = Decoration_class([],[],[],[]);

            List = Decoration.getListOfLines(self.msg, self.newLine);

            lenList = length(List);
            for i = 1 : lenList
                ListJustified       = Decoration.wrapText(List{i}, self.wrapWidth, self.wrapSplit , []);
                lenListJustified    = length(ListJustified);
                for ijustified = 1 : lenListJustified
                    padTopCurrent = 0;
                    padBotCurrent = 0;
                    if (i == 1) && (ijustified == 1)                        , padTopCurrent = self.marginTop; end % the very first line
                    if (i == lenList) && (ijustified == lenListJustified)   , padBotCurrent = self.marginBot; end % the very last line
                    Decoration.write(self.outputUnit, padTopCurrent, padBotCurrent, 1, self.fullprefix + ListJustified{ijustified});
                end
            end

            if self.marginBot==2, Decoration.write(self.outputUnit,[],[],[],[]); end

            if self.resetEnabled
                self.reset();
            end

        end % function informUser

        %***************************************************************************************************************************
        %***************************************************************************************************************************

        function reset(self)

            self.msg        = "";
            self.prefix     = "";
            self.wrapSplit  = " ";
            self.wrapWidth  = 100;
            self.marginTop  = 1;
            self.marginBot  = 1;
            self.outputUnit = 1;
            self.occurred   = false;
            self.newLine    = newline;
            self.stat       = -intmax;
            self.statNull   = -intmax;

        end

        %***************************************************************************************************************************
        %***************************************************************************************************************************

    end % methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    methods (Hidden)
        % These methods have been implemented to override the default 'handle' class methods,
        % so that they won't pop-up after pressing 'Tab' button.
        function addlistener    ();  end
        function delete         ();  end
        function findobj        ();  end
        function findprop       ();  end
        function valid          ();  end
        function listener       ();  end
        function notify         ();  end
    end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

end