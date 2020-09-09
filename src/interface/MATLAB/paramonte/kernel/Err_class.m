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

classdef Err_class < handle

    properties (Constant)
        CLASS_NAME              = "@Err_class"
        ERR_HANDLING_REQUESTED  = false
    end

    properties
        msg
        prefix
        advance
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
            self.advance = true;
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
                Decoration.write(1, 0, 0, 1, self.fullprefix + "For more information, please see the output report file.");
                Decoration.write(1, 0, 2, 1, self.fullprefix + "Gracefully Exiting on image " + strtrim(imageChar) + ".");
            else
                Decoration.write(self.outputUnit, 1, 0, 1, self.fullprefix + "Please Correct the error(s) and rerun the simulation.");
                Decoration.write(self.outputUnit, 1, 0, 1, self.fullprefix + "For further help, contact the ParaMonte authors via:");
                Decoration.write(self.outputUnit, 0, 0, 1, self.fullprefix + "cdslab.org/ParaMonte/");
                Decoration.write(self.outputUnit, 1, 2, 1, self.fullprefix + "Gracefully Exiting on image " + strtrim(imageChar) + ".");
            end
            error(newline+"Fatal error occurred."+newline);
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
            Decoration.advance = self.advance;

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