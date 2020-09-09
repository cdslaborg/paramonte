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

classdef Decoration_class < handle

    properties (Constant)
        CLASS_NAME                  = '@Decoration_mod'
        DECORATION_WIDTH            = 132
        DECORATION_THICKNESS_HORZ   = 4
        DECORATION_THICKNESS_VERT   = 1
        STAR                        = "*"
        TAB                         = "    "
    end

    properties
        tab                         = []
        text                        = []
        symbol                      = []
        advance                     = []
        ListOfLines                 = {}
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = Decoration_class(tabStr, symbol, text, ListOfLines)
            if isempty(tabStr), self.tab = self.TAB;        else, self.tab = tabStr;    end
            if isempty(symbol), self.symbol = self.STAR;    else, self.symbol = symbol; end
            self.advance = true;
            self.text = text;
            self.ListOfLines = ListOfLines;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function writeDecoratedText(self, text, symbol, width, thicknessHorz, thicknessVert, marginTop, marginBot, outputUnit, newLine)
            if isempty(thicknessVert)
                thicknessVertDefault = self.DECORATION_THICKNESS_VERT;
            else
                thicknessVertDefault = thicknessVert;
            end

            if isempty(newLine)
                self.write(outputUnit,   marginTop,    0,          thicknessVertDefault,   self.drawLine(symbol,width)                      );
                self.write(outputUnit,   0,            0,          1,                      self.sandwich(text,symbol,width,thicknessHorz)   );
                self.write(outputUnit,   0,            marginBot,  thicknessVertDefault,   self.drawLine(symbol,width)                      );
            else
                self.writeDecoratedList(self.getListOfLines(text,newline), symbol, width, thicknessHorz, thicknessVert, marginTop, marginBot, outputUnit);
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function writeDecoratedList(self, ListOfLines, symbol, width, thicknessHorz, thicknessVert, marginTop, marginBot, outputUnit)

            if isempty(thicknessVert)
                thicknessVertDefault = self.DECORATION_THICKNESS_VERT;
            else
                thicknessVertDefault = thicknessVert;
            end

            self.write( outputUnit, marginTop, 0, thicknessVertDefault, self.drawLine(symbol,width) );

            for i = 1 : length(ListOfLines)
                self.write( outputUnit, 0, 0, 1, self.sandwich(ListOfLines{i}, symbol, width, thicknessHorz) );
            end

            self.write( outputUnit, 0, marginBot, thicknessVertDefault, self.drawLine(symbol,width) );

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function line = drawLine(self, symbol, width)

            symbol = convertStringsToChars(symbol);

            if isempty(symbol)
                decorationSymbol = self.STAR;
            else
                if length(symbol) < 1
                decorationSymbol = " ";
                else
                decorationSymbol = symbol;
                end
            end

            if isempty(width)
                decorationWidth = self.DECORATION_WIDTH;
            else
                decorationWidth = width;
            end

            symbolLen = length(decorationSymbol);
            decorationSymbol = convertStringsToChars(decorationSymbol);
            line = decorationSymbol;

            for i = 1 : floor(decorationWidth / symbolLen)-1
                line = strcat(line, decorationSymbol);
            end
            line = strcat(line, decorationSymbol(1:mod(decorationWidth, length(decorationSymbol))));

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function sandwichedText = sandwich(self, text, symbol, width, thicknessHorz)
            text    = convertStringsToChars(text);
            symbol  = convertStringsToChars(symbol);

            if isempty(symbol),         decorationSymbol = self.STAR;                               else, decorationSymbol = symbol;                end
            if isempty(width),          decorationWidth = self.DECORATION_WIDTH;                    else, decorationWidth = width;                  end
            if isempty(thicknessHorz),  decorationThicknessHorz = self.DECORATION_THICKNESS_HORZ;   else, decorationThicknessHorz = thicknessHorz;  end
            if isempty(text),           decorationText = "";                                        else, decorationText = text;                    end

            decorationText      = convertStringsToChars(decorationText);
            decorationSymbol    = convertStringsToChars(decorationSymbol);

            decorationTextLen   = length(decorationText);
            decorationSymbolLen = length(decorationSymbol);
            sandwichedText      = blanks(decorationWidth);

            if decorationWidth < 1, sandwichedText = ""; return; end
            %-----------------------------------------------------------------------------------------------------------------------
            % decorationWidth setup
            if decorationSymbolLen >= decorationThicknessHorz
                for i = 1 : decorationThicknessHorz
                    sandwichedText(i) = decorationSymbol(i);
                end
                j = 1;
                for i = decorationWidth - decorationThicknessHorz + 1 : decorationWidth
                    sandwichedText(i) = decorationSymbol(j);
                    j = j + 1;
                end
            else
                j = 1;
                for i = 1 : decorationThicknessHorz
                    sandwichedText(i) = decorationSymbol(j);
                    if j == decorationSymbolLen, j = 0; end
                    j = j + 1;
                end
                j = 1;
                for i = decorationWidth - decorationThicknessHorz + 1 : decorationWidth
                    sandwichedText(i) = decorationSymbol(j);
                    if j == decorationSymbolLen, j = 0; end
                    j = j + 1;
                end
            end
            %-----------------------------------------------------------------------------------------------------------------------
            % sandwichedText setup
            if decorationTextLen > (decorationWidth - 2*decorationThicknessHorz)
                for i = decorationThicknessHorz + 1 : decorationWidth - decorationThicknessHorz
                    sandwichedText(i) = decorationText(i-decorationThicknessHorz);
                end
            else
                j = 1;
                for i = round(decorationWidth/2) - round(decorationTextLen/2) + 1 : round(decorationWidth/2) - round(decorationTextLen/2) + decorationTextLen
                    sandwichedText(i) = decorationText(j);
                    j = j + 1;
                end
            end
        end % function sandwichedText

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function write(self, outputUnit, marginTop, marginBot, count, string)
            if self.advance; stringFormat = "%s\n"; else stringFormat = "%s"; end
            if ~isempty(outputUnit), logFileUnit = outputUnit; else, logFileUnit = 1; end
            %-----------------------------------------------------------------------------------------------------------------------
            if ~isempty(marginTop)
                for i = 1 : marginTop
                    fprintf(logFileUnit, "\n");
                end
            end
            %-----------------------------------------------------------------------------------------------------------------------
            if ~isempty(count)
                thisManyTimes = count;
            else
                thisManyTimes = 1;
            end
            %-----------------------------------------------------------------------------------------------------------------------
            if ~isempty(string)
                for i = 1 : thisManyTimes
                    fprintf(logFileUnit, stringFormat, string);
                end
            elseif ~(~isempty(marginBot) && ~isempty(marginTop))
                for i = 1 : thisManyTimes
                    fprintf(logFileUnit, "\n");
                end
            end
            %-----------------------------------------------------------------------------------------------------------------------
            if ~isempty(marginBot)
                for i = 1 : marginBot
                    fprintf(logFileUnit, "\n");
                end
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function ListOfLines = wrapText(self, string, width, delim, pad)
            FUNCTION_NAME = [self.CLASS_NAME,'@wrapText()'];

            %string = strrep(string, newline + newline, newline + ' ' + newline);

            if isempty(width), width = 132; end

            pad     = convertStringsToChars(pad);
            delim   = convertStringsToChars(delim);
            string  = convertStringsToChars(string);

            padLen      = length(pad);
            delimLen    = length(delim);
            stringLen   = length(string);

            if stringLen == 0
                ListOfLines{1} = '';
                return
            elseif delimLen >= stringLen
                ListOfLines{1} = string;
                return
            elseif delimLen == 0 % enforce wrapping at any word as necessary
                lineCount = floor(stringLen / width) + 1;
                ListOfLines = cell(lineCount,1);
                for indx = 1 : lineCount - 1;
                    ListOfLines{indx} = string(width*(indx-1)+1 : width*indx);
                end
                ListOfLines{lineCount} = string(width*(lineCount-1)+1:stringLen);
                return
            end

            % get the initial pad size, and the locations of delim ends.

            IsEndOfDelimLoc = zeros(stringLen,1);
            IsEndOfDelimLoc = 0;
            istart = 1;
            iend = istart + delimLen - 1;
            padLength = 0;
            isPadZone = true;
            if  padLen==0
                isPadZone = false;
            end
            while true % blockFindDelim
                if  iend==stringLen
                    IsEndOfDelimLoc(stringLen) = 1;
                    break % blockFindDelim
                end
                if  isPadZone && mod(iend,padLen)==0 && string(istart:iend)==pad
                    padLength = iend;
                else
                    isPadZone = false;
                end
                if  string(istart:iend)==delim
                    IsEndOfDelimLoc(iend) = 1;
                else
                    IsEndOfDelimLoc(iend) = 0;
                end
                istart = istart + 1;
                iend = iend + 1;
            end % blockFindDelim

            % create a vector of delim-end indices

            numDelimEndLoc = sum(IsEndOfDelimLoc);
            if  numDelimEndLoc==0
                ListOfLines{1} = string;
                return
            else
                EndOfDelimLoc = zeros(numDelimEndLoc,1);
                counter = 0;
                for indx = 1:stringLen
                    if  IsEndOfDelimLoc(indx)==1
                        counter = counter + 1;
                        EndOfDelimLoc(counter) = indx;
                    end
                end
            end
            EndOfDelimLoc = EndOfDelimLoc(1:counter);

            % compute the number wrappings to be done

            EndOfLineLoc = zeros(numDelimEndLoc+2,1); % (0:numDelimEndLoc+1) ) % consider the maximum possible number of lines
            lineCount = 0;
            lineCountPlusOne = lineCount + 1;
            padLengthDynamic = 0; % first wrap does not need padding
            indxOld = 1;
            indx = 0;
            oldLineLen = -intmax;
            while true % blockFindLine
                indx = indx + 1;
                if  indx>numDelimEndLoc
                    break % blockFindLine
                end
                newLineLen = padLengthDynamic + EndOfDelimLoc(indx)-EndOfLineLoc(lineCountPlusOne);
                if  newLineLen<=width
                    oldLineLen = newLineLen;
                    continue % blockFindLine
                else
                    % swap the commented block with the uncommented to switch from better to good wrapping style.
                    lineCount = lineCount + 1;
                    lineCountPlusOne = lineCountPlusOne + 1;
                    if  indx-1>indxOld; % ensure there is at least one delim before the wrapping point
                        % comment the following line to keep the max line length, strictly less than width (if possible).
                        if width-oldLineLen > newLineLen-width
                            indx = indx + 1; % removing the last token would make the line more beautiful
                        end
                        EndOfLineLoc(lineCountPlusOne) = EndOfDelimLoc(indx-1);
                    else
                        EndOfLineLoc(lineCountPlusOne) = EndOfDelimLoc(indx);
                    end
                    indxOld = indx;
                    padLengthDynamic = padLength;
                end
            end % blockFindLine

            % add the remaining end of the string as a separate line

            if  EndOfLineLoc(lineCountPlusOne) < stringLen || lineCount == 0
                lineCount = lineCount + 1;
                lineCountPlusOne = lineCountPlusOne + 1;
                EndOfLineLoc(lineCountPlusOne) = stringLen;
            end
            EndOfLineLoc = EndOfLineLoc(EndOfLineLoc~=0);

            % ensure the line count makes sense

            if  lineCount ~= length(EndOfLineLoc)
                error   (   [ FUNCTION_NAME ...
                            , ': Fatal error occurred. lineCount ~= length(EndOfLineLoc): ' ...
                            , num2str(lineCount), '~=', num2str(length(EndOfLineLoc)) ...
                            ] ...
                        );
            end

            % construct the wrappings

            ListOfLines = cell(lineCount,1);
            indx = 1;
            ListOfLines{indx} = string(1:EndOfLineLoc(indx));
            for indx = 2:lineCount
                if  padLength==0 && EndOfLineLoc(indx-1)+1>EndOfLineLoc(indx)
                    disp(string)
                end
                ListOfLines{indx} = [ string(1:padLength), string(EndOfLineLoc(indx-1)+1:EndOfLineLoc(indx)) ];
            end

        end % function wrapText

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods(Static)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function ListOfLines = getListOfLines(string, delimiter)

            if isempty(delimiter)
                ListOfLines{1} = string;
                return
            end

            string      = convertStringsToChars(string);
            delimiter   = convertStringsToChars(delimiter);

            stringLen   = length(string);
            delimLen    = length(delimiter);

            delimLenMinusOne = delimLen - 1;

            if (delimLen == 0) || (stringLen == 0) || (stringLen < delimLen)
                ListOfLines{1} = string;
                return
            end

            delimIsCStyle   = ((delimLen == 2) && (delimiter=="\n"));

            maxNumSplit     = 1 + floor(stringLen / delimLen);
            ListOfLines     = {blanks(maxNumSplit)};

            counterLine     = 0;
            counterRecord   = 0;
            counterString   = 1;

            while true
                if (counterString + delimLenMinusOne) > stringLen
                    counterLine = counterLine + 1;
                    if counterRecord == 0
                        ListOfLines{counterLine} = string(counterString : stringLen);
                    else
                        ListOfLines{counterLine} = dumstr(1 : counterRecord) + convertCharsToStrings(string(counterString : stringLen));
                    end
                    break;
                end
                if string(counterString : counterString + delimLenMinusOne) == delimiter
                    counterLine = counterLine + 1;
                    if counterRecord == 0
                        ListOfLines{counterLine} = "";
                    else
                        ListOfLines{counterLine} = dumstr(1 : counterRecord);
                        counterRecord = 0;
                    end
                    counterString = counterString + delimLen;
                    if counterString > stringLen
                        counterLine = counterLine + 1;
                        ListOfLines{counterLine} = "";
                        break;
                    end
                elseif delimIsCStyle && (string(counterString : counterString) == '\')
                    counterString = counterString + 1;
                    counterRecord = counterRecord + 1;
                    dumstr(counterRecord : counterRecord) = '\';
                    if string(counterString : counterString) == '\', counterString = counterString + 1; end
                else
                    counterRecord = counterRecord + 1;
                    dumstr(counterRecord : counterRecord) = string(counterString : counterString);
                    counterString = counterString + 1;
                end
                ListOfLines = ListOfLines(1 : counterLine);
            end

        end % function ListOfLines

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end


%***********************************************************************************************************************************
%***********************************************************************************************************************************

end