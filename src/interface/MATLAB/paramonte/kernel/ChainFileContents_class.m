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

classdef ChainFileContents_class < handle

    properties (Constant)
        CLASS_NAME          = "@ParaDRAMChainFileContents_mod"
        NUM_DEF_COL         = 7 % number of columns in the chain file other than the State columns
        COL_HEADER_DEFAULT  =   [ "ProcessID            "   ...
                                , "DelayedRejectionStage"   ...
                                , "MeanAcceptanceRate   "   ...
                                , "AdaptationMeasure    "   ...
                                , "BurninLocation       "   ...
                                , "SampleWeight         "   ...
                                , "SampleLogFunc        "   ...
                                ]
                                                                  
    end

    properties
        ndim        = 0
        lenHeader   = 0
        numDefCol   = ChainFileContents_class.NUM_DEF_COL
        Count       = Count_class()
        ProcessID   = []                    % the vector of the ID of the images whose function calls haven been accepted
        DelRejStage = []                    % delayed rejection stages at which the proposed states were accepted
        Adaptation  = []                    % the vector of the adaptation measures at the MCMC accepted states.
        MeanAccRate = []                    % the vector of the average acceptance rates at the given point in the chain
        BurninLoc   = []                    % the burnin locations at the given locations in the chains
        Weight      = []                    % the vector of the weights of the MCMC accepted states.
        LogFunc     = []                    % the vector of LogFunc values corresponding to the MCMC states.
        State       = []                    % the (nd,chainSize) MCMC chain of accepted proposed states
        ColHeader   = [""]                  % column headers of the chain file
        delimiter   = []                    % delimiter used to separate objects in the chain file
        Err         = Err_class()

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        % if chainFilePath is given, the rest of the optional arguments must be also given
        function self = ChainFileContents_class(ndim, variableNameList, chainFilePath, chainSize, chainFileForm, lenHeader, delimiter, targetChainSize)

            Err             = Err_class();
            Err.occurred    = false;
            self.ndim       = ndim;

            % set up the chain file column header

            for icol = 1 : ChainFileContents_class.NUM_DEF_COL
                self.ColHeader(icol) = strtrim(self.COL_HEADER_DEFAULT(icol));
            end

            if ~isempty(variableNameList)
                for icol = ChainFileContents_class.NUM_DEF_COL + 1 : ChainFileContents_class.NUM_DEF_COL + ndim
                    self.ColHeader(icol) = strtrim(variableNameList{icol - ChainFileContents_class.NUM_DEF_COL});
                end
            end

            % set up other variables if given

            if ~isempty(lenHeader),         self.lenHeader      = lenHeader;        end
            if ~isempty(delimiter),         self.delimiter      = delimiter;        end
            if ~isempty(targetChainSize),  self.Count.target  = targetChainSize; end

            % read the chain file if the path is given

            if ~isempty(chainFilePath), Err = self.get(chainFilePath, chainFileForm, chainSize, lenHeader, ndim, delimiter, targetChainSize); end
            if Err.occurred
                self.Err.occurred   = true;
                self.Err.msg        = Err.msg;
                return
            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = get(self, chainFilePath, chainFileForm, chainSize, lenHeader, ndim, delimiter, targetChainSize)
            % returns the contents of a ParaDRAM simulation output chain file, always in compact format, regardless of chainFileFormat.
            % targetChainSize: >= chainSize, used for the allocation of the chain components
            % chainSize: <= targetChainSize, the first chainSize elements of the Chain components will contain the Chain information read from the chain file.
            %                               , the Chain component elements beyond chainSize will be set to zero.

            FUNCTION_NAME = self.CLASS_NAME + "@get()";

            Err                     = Err_class();
            Err.occurred            = false;
            chainFilePathTrimmed    = strtrim(chainFilePath);

            % inquire(file=chainFilePathTrimmed,exist=fileExists,opened=fileIsOpen,number=chainFileUnit)

            fileExists              = isfile(chainFilePathTrimmed);
            chainFileUnit           = fopen(chainFilePathTrimmed);

            if fileExists   % blockFileExistence

                % set up chain file format

                isBinary    = false;
                isCompact   = false;
                isVerbose   = false;

                if lower(chainFileForm)     == "binary"
                    isBinary    = true;
                elseif lower(chainFileForm) == "compact"
                    isCompact   = true;
                elseif lower(chainFileForm) == "verbose"
                    isVerbose   = false;
                else
                    Err.occurred    = true;
                    Err.msg         = FUNCTION_NAME + ": Unrecognized chain file form: " + chainFileForm;
                    return
                end

                if isBinary
                    thisForm = "unformatted";
                    if isempty(ndim) || isempty(lenHeader) || isempty(delimiter)
                        Err.occurred    = true;
                        Err.msg         = FUNCTION_NAME + ": If the chain file is in binary form, chainSize, lenHeader, delimiter, and ndim must be provided by the user.";
                        return
                    end
                else
                    thisForm = "formatted";
                end

                if chainFileUnit == -1
                    Err.occurred    = true;
                    Err.msg         = FUNCTION_NAME + ": Unable to open the file located at: " + chainFilePathTrimmed + newline;
                    return
                else
%                   fclose(chainFileUnit);  % xxx
                end

                % get the number of records in file, minus header line

%                if ~isempty(chainSize)
%                    chainSizeDefault = chainSize;
%                else   % here chainSizeDefault is indeed max(chainSize) depending on the file format: verbose or compact
                    if isBinary
                        [chainFileUnit, Err.msg] = fopen(chainFilePathTrimmed);
                        if Err.msg
                            Err.occurred    = true;
                            Err.msg         = FUNCTION_NAME + ": Unable to open the file located at: " + chainFilePathTrimmed + newline;
                            return
                        end

                        Record.value = fread(chainFileUnit);

                        State(ndim) = zeros(ndim);

                        chainSizeDefault = 0;

                        while true  % loopFindChainSizeDefault
                            % read(chainFileUnit,iostat=Err%stat) processID, delRejStage, meanAccRate, adaptation, burninLoc, weight, logFunc, State
                            if Err.stat == 0
                                chainSizeDefault = chainSizeDefault + 1;
                            elseif is_iostat_end(Err.stat)
                                break % loopFindChainSizeDefault
                            elseif is_iostat_eor(Err.stat)
                                Err.occurred    = true;
                                Err.msg         = FUNCTION_NAME + ": Incomplete record detected while reading the input binary chain file at: " + chainFilePathTrimmed + newline;
                                return
                            else
                                Err.occurred    = true;
                                Err.msg         = FUNCTION_NAME + ": IO error occurred while reading the input binary chain file at: " + chainFilePathTrimmed + newline;
                                return
                            end
                        end
                        fclose(chainFileUnit);
                    else    % is not binary
                        chainSizeDefault = getNumRecordInFile(chainFileUnit);
                        if chainSizeDefault == 0, self.Count.compact = 0; return; end
                        chainSizeDefault = chainSizeDefault - 1;    % subtract header
                    end
 %               end

                % set the number of elements in the Chain components

                if ~isempty(targetChainSize)    % in restart mode, this must be always the case
                    self.Count.target = targetChainSize;
                else
                    self.Count.target = chainSizeDefault;
                end
                % if self.Count.target < chainSizeDefault
                    % Err.occurred    = true;
                    % Err.msg         = FUNCTION_NAME + ": input targetChainSize cannot be smaller than the input chainSize:" + newline  ...
                                    % + "    targetChainSize = " + num2str(self.Count.target) + newline                                 ...
                                    % + "           chainSize = " + num2str(chainSizeDefault) + newline                                   ...
                                    % + "It appears that the user has manipulated the output chain file."                                 ...
                                    % ;
                    % return
                % end

                % allocate Chain components

                self.ProcessID      = zeros(1, self.Count.target);
                self.DelRejStage    = zeros(1, self.Count.target);
                self.MeanAccRate    = zeros(1, self.Count.target);
                self.Adaptation     = zeros(1, self.Count.target);
                self.BurninLoc      = zeros(1, self.Count.target);
                self.Weight         = zeros(1, self.Count.target);
                self.LogFunc        = zeros(1, self.Count.target);

                % find the delimiter

                if ~isempty(delimiter)  % blockFindDelim

                    self.delimiter = delimiter;

                else % blockFindDelim

                    [chainFileUnit, Err.msg] = fopen(chainFilePathTrimmed);
                    if Err.msg
                        Err.occurred    = true;
                        Err.msg         = FUNCTION_NAME + ": Unable to open the file located at: " + chainFilePathTrimmed + "." + newline;
                        return
                    end

                    fgets(chainFileUnit);   % skip the header
                    Record.value = fgets(chainFileUnit);    % read the first numeric row in string format
                    fclose(chainFileUnit);

                    Record.value    = strtrim(Record.value);
                    delimHasEnded   = false;
                    delimHasBegun   = false;
                    delimiterLen    = 0;

                    for i = 1 : length(Record.value) - 1  % loopSearchDelimiter
                        if Record.isDigit(Record.value(i:i))
                            if delimHasBegun, delimHasEnded = true; end
                        elseif Record.value(i:i) == "." || Record.value(i:i) == "+" || Record.value(i:i) == "-"
                            if delimHasBegun
                                delimHasEnded = true;
                            else
                                Err.occurred    = true;
                                Err.msg         = FUNCTION_NAME + ": The file located at: " + chainFilePathTrimmed + newline                                                    ...
                                                + "has unrecognizable format. Found " + Record.value(i:i) + " in the first column, while expecting positive integer." + newline ...
                                                ;
                                return
                            end
                        else
                            if i==1 % here it is assumed that the first column in chain file always contains integers
                                Err.occurred    = true
                                Err.msg         = FUNCTION_NAME + ": The file located at: " + chainFilePathTrimmed + newline    ...
                                                + "has unrecognizable format." + newline                                        ...
                                                ;
                                return
                            else
                                delimHasBegun   = true;
                                delimiterLen    = delimiterLen + 1;
                                self.delimiter(delimiterLen:delimiterLen) = Record.value(i:i);
                            end
                        end
                        if delimHasEnded, break; end
                    end %loopSearchDelimiter

                    if ~(delimHasBegun && delimHasEnded)
                        Err.occurred    = true;
                        Err.msg         = FUNCTION_NAME + ": The file located at: " + chainFilePathTrimmed + newline + "has unrecognizable format. Could not identify the column delimiter." + newline;
                        return
                    else
                        self.delimiter = strtrim(self.delimiter(1:delimiterLen));
                        delimiterLen = length(self.delimiter);
                        if delimiterLen == 0
                            self.delimiter  = " ";
                            delimiterLen    = 1;
                        end
                    end

                end % blockFindDelim

                % find the number of dimensions of the state (the number of function variables)

                if ~isempty(ndim)
                    self.ndim = ndim;
                else
                    Record.Parts = Record.SplitStr(Record.value, self.delimiter, Record.nPart); %xxx
                    self.numDefCol = 0;
                    for i = 1 : Record.nPart    % loopFindNumDefCol
                        if index(Record.Parts(i), "LogFunc") > 0
                            self.numDefCol = i;
                            break
                        end
                    end % loopFindNumDefCol
                    if self.numDefCol ~= self.NUM_DEF_COL || self.numDefCol == 0
                        Err.occurred    = true;
                        Err.msg         = FUNCTION_NAME + ": Internal error occurred. CFC.numDefCol ~= CFC.NUM_DEF_COL: " + num2str(self.numDefCol) + num2str(self.NUM_DEF_COL);
                        return
                    end
                    self.ndim = Record.nPart - self.NUM_DEF_COL;
                end

                % reopen the file to read the contents

                [chainFileUnit, Err.msg] = fopen(chainFilePathTrimmed);
                if Err.msg
                    Err.occurred    = true;
                    Err.msg         = FUNCTION_NAME + ": Unable to open the file located at: " + chainFilePathTrimmed  + "." + newline;
                    return
                end

                % first read the column headers

                if isBinary
                  % Record.value = fgets(chainFileUnit);
                    Record.value = convertCharsToStrings(Record.value);
                else
                    Record.value = fgets(chainFileUnit);
                    Record.value = convertCharsToStrings(Record.value);
                end
                
                self.ColHeader = split(strtrim(Record.value), self.delimiter);
                for i = 1 : length(self.ColHeader)
                    self.ColHeader(i) = strtrim(self.ColHeader(i));
                end

                % read the chain

                self.State          = zeros(self.ndim, self.Count.target);
                self.Count.verbose  = 0;

                if isBinary
                    % do iState = 1, chainSizeDefault
                        % read(chainFileUnit  ) CFC%ProcessID                (iState)    &
                                            % , CFC%DelRejStage              (iState)    &
                                            % , CFC%MeanAccRate              (iState)    &
                                            % , CFC%Adaptation               (iState)    &
                                            % , CFC%BurninLoc                (iState)    &
                                            % , CFC%Weight                   (iState)    &
                                            % , CFC%LogFunc                  (iState)    &
                                            % , CFC%State         (1:CFC%ndim,iState)
                        % CFC%Count%verbose = CFC%Count%verbose + CFC%Weight(iState)
                    % end do
                elseif isCompact
                    for iState = 1 : chainSizeDefault   % xxx
                        Record.value = fgets(chainFileUnit);
                        if Record.value ~= -1   % xxx
                            Record.Parts = split(strtrim(Record.value), self.delimiter);
                            Record.nPart = length(Record.Parts);
                            self.ProcessID  (iState) = str2num(Record.Parts{1});
                            self.DelRejStage(iState) = str2num(Record.Parts{2});
                            self.MeanAccRate(iState) = str2num(Record.Parts{3});
                            self.Adaptation (iState) = str2num(Record.Parts{4});
                            self.BurninLoc  (iState) = str2num(Record.Parts{5});
                            self.Weight     (iState) = str2num(Record.Parts{6});
                            self.LogFunc    (iState) = str2num(Record.Parts{7});
                            for i = 1 : self.ndim
                                self.State(i, iState) = str2num(Record.Parts{7+i});
                            end
                        end
                        self.Count.verbose = self.Count.verbose + self.Weight(iState);
                    end
                else % is verbose form
                    if chainSizeDefault > 0 % blockChainSizeDefault
                        self.Count.compact = 1;
                        State = zeros(1, ndim);
                        % read the first sample
                        Record.value = fgets(chainFileUnit);
                        Record.Parts = split(strtrim(Record.value), self.delimiter);
                        self.ProcessID  (self.Count.compact) = str2num(Record.Parts{1});
                        self.DelRejStage(self.Count.compact) = str2num(Record.Parts{2});
                        self.MeanAccRate(self.Count.compact) = str2num(Record.Parts{3});
                        self.Adaptation (self.Count.compact) = str2num(Record.Parts{4});
                        self.BurninLoc  (self.Count.compact) = str2num(Record.Parts{5});
                        self.Weight     (self.Count.compact) = str2num(Record.Parts{6});
                        self.LogFunc    (self.Count.compact) = str2num(Record.Parts{7});
                        for i = 1 : self.ndim
                            self.State(i, self.Count.compact) = str2num(Record.Parts{7+i});
                        end

                        % read the rest of samples beyond the first, if any exist
                        newUniqueSampleDetected = false;
                        for iState = 2 : chainSizeDefault   % loopOverChainfFileContents
                            Record.value = fgets(chainFileUnit);
                            if Record.value ~= -1
                                Record.Parts = split(strtrim(Record.value), self.delimiter);
                                ProcessID   = str2num(Record.Parts{1});
                                DelRejStage = str2num(Record.Parts{2});
                                MeanAccRate = str2num(Record.Parts{3});
                                Adaptation  = str2num(Record.Parts{4});
                                BurninLoc   = str2num(Record.Parts{5});
                                Weight      = str2num(Record.Parts{6});
                                LogFunc     = str2num(Record.Parts{7});
                                for i = 1 : self.ndim
                                    State(i) = str2num(Record.Parts{7+i});
                                end
                            end
                            % increment CFC%Count%compact if new sample detected
                            newUniqueSampleDetected =  LogFunc      ~= self.LogFunc    (self.Count.compact);
                                                % || MeanAccRate    ~= self.MeanAccRate(self.Count.compact)   ...
                                                % || Adaptation     ~= self.Adaptation (self.Count.compact)   ...
                                                % || BurninLoc      ~= self.BurninLoc  (self.Count.compact)   ...
                                                % || Weight         ~= self.Weight     (self.Count.compact)   ...
                                                % || DelRejStage    ~= self.DelRejStage(self.Count.compact)   ...
                                                % || ProcessID      ~= self.ProcessID  (self.Count.compact)   ...
                                                % || any( self.State(1:self.ndim, self.Count.compact) ~= self.State(1:self.ndim, self.Count.compact) );
                            if newUniqueSampleDetected
                                self.Count.compact  =  self.Count.compact + 1;
                                if self.Count.target < self.Count.compact
                                    Err.occurred    = true;
                                    Err.msg         = FUNCTION_NAME + ": Fatal error occurred. CFC.Count.target < CFC.Count.compact: "          ...
                                                    + num2str(self.Count.target) + " ~= " + num2str(self.Count.compact)                         ...
                                                    + "The contents of the input chain file is longer than the user-requested allocation size." ...
                                                    ;
                                    return
                                end
                            else
                                Weight = self.Weight(self.Count.compact) + 1;
                            end
                            
                            % write the latest sample
                            self.LogFunc         (self.Count.compact)   = LogFunc           ;
                            self.MeanAccRate     (self.Count.compact)   = MeanAccRate       ;
                            self.Adaptation      (self.Count.compact)   = Adaptation        ;
                            self.BurninLoc       (self.Count.compact)   = BurninLoc         ;
                            self.Weight          (self.Count.compact)   = Weight            ;
                            self.DelRejStage     (self.Count.compact)   = DelRejStage       ;
                            self.ProcessID       (self.Count.compact)   = ProcessID         ;
                            self.State(1:self.ndim, self.Count.compact) = State(1:self.ndim);
                        end % loopOverChainfFileContents
                    else    % blockChainSizeDefault
                        self.Count.compact = 0;
                        self.Count.verbose = 0;
                    end
                end

                if isBinary || isCompact
                    self.Count.compact  = chainSizeDefault;
                else
                    self.Count.verbose  = chainSizeDefault;
                    if self.Count.verbose ~= sum(self.Weight(1:self.Count.compact))
                        Err.occurred    = true;
                        Err.msg         = FUNCTION_NAME + ": Internal error occurred. CountVerbose ~= sum(Weight): "            ...
                                        + num2str(self.Count.verbose) + " ~= "+ num2str(sum(self.Weight(1:self.Count.compact))) ...
                                        ;
                        return
                    elseif isempty(targetChainSize)
                        self.ProcessID      = self.ProcessID   (1:self.Count.compact);
                        self.DelRejStage    = self.DelRejStage (1:self.Count.compact);
                        self.MeanAccRate    = self.MeanAccRate (1:self.Count.compact);
                        self.Adaptation     = self.Adaptation  (1:self.Count.compact);
                        self.BurninLoc      = self.BurninLoc   (1:self.Count.compact);
                        self.Weight         = self.Weight      (1:self.Count.compact);
                        self.LogFunc        = self.LogFunc     (1:self.Count.compact);
                        self.State          = self.State       (1:self.ndim, 1:self.Count.compact);
                    end
                end

                fclose(chainFileUnit);

                % set the rest of elements to null values
                if self.Count.target > chainSizeDefault
                    % self.nullify(startIndex = self.Count.compact + 1, endIndex = self.Count.target);
                end

            else % blockFileExistence

                Err.occurred    = true;
                Err.msg         = FUNCTION_NAME + ": The chain file does not exist in the given file path: " + chainFilePathTrimmed;
                return

            end % blockFileExistence

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function getLenHeader(CFC, ndim, isBinary, chainFileFormat)

            FUNCTION_NAME = self.CLASS_NAME + "@getLenHeader()";

            self.Err.occurred = false;
            if isBinary
              % fprintf( record , "(*(g0,:,','))" ) (CFC%ColHeader(i)%record, i=1,CFC%numDefCol+ndim)
                record = "";
                for i = 1 : self.numDefCol+ndim, record = record + self.ColHeader(i); end
            else
                if ~isempty(chainFileFormat)
                  % fprintf( record,  chainFileFormat,) (CFC%ColHeader(i)%record, i=1,CFC%numDefCol+ndim)
                    for i = 1 : self.numDefCol+ndim, record = record + self.ColHeader(i); end
                else
                    self.Err.occurred   = true;
                    self.Err.msg        = FUNCTION_NAME + "Internal error occurred. For formatted chain files, chainFileFormat must be given.";
                    self.Err.abort();
                end
            end

            self.lenHeader = length(strtrim(record));

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function writeHeader(self, ndim, chainFileUnit, isBinary, headerFormat)

            FUNCTION_NAME   = self.CLASS_NAME + "@writeHeader()";

            self.Err.occurred = false;
            if isBinary
                % write( record , "(*(g0,:,','))" ) (CFC%ColHeader(i)%record, i=1,CFC%numDefCol+ndim)
                % write(chainFileUnit) trim(adjustl(record))
                % deallocate(record)
            else
                fprintf(chainFileUnit, headerFormat, self.ColHeader);
            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function writeChainFile(self, ndim, compactStartIndex, compactEndIndex, chainFileUnit, chainFileForm, chainFileFormat, headerFormat)

            FUNCTION_NAME       = self.CLASS_NAME + "@writeChainFile()";

            self.Err.occurred   = false;

            isBinary            = false;
            isVerbose           = false;
            isCompact           = false;

            if     chainFileForm == "binary"
                isBinary    = true;
            elseif chainFileForm == "compact"
                isCompact   = true;
            elseif chainFileForm == "verbose"
                isVerbose   = true;
            else
                self.Err.occurred   = true;
                self.Err.msg        = FUNCTION_NAME + "Internal error occurred. Unknown chain file format: " + chainFileForm;
            end

            if ~isBinary && isempty(chainFileFormat)
                self.Err.occurred   = true;
                self.Err.msg        = FUNCTION_NAME + "Internal error occurred. For formatted chain files, chainFileFormat must be given."
            end
            if self.Err.occurred, self.Err.abort([], [], []); end

            self.writeHeader(ndim, chainFileUnit, isBinary, headerFormat);

            if compactStartIndex <= compactEndIndex
                if isCompact
                    for i = compactStartIndex : compactEndIndex
                        fprintf(chainFileUnit   , chainFileFormat                   ...
                                                    , self.ProcessID    (i)         ...
                                                    , self.DelRejStage  (i)         ...
                                                    , self.MeanAccRate  (i)         ...
                                                    , self.Adaptation   (i)         ...
                                                    , self.BurninLoc    (i)         ...
                                                    , self.Weight       (i)         ...
                                                    , self.LogFunc      (i)         ...
                                                    , self.State (1:ndim,i)         ...
                                                    ) ;
                    end
                elseif isBinary
                    % for i = compactStartIndex : compactEndIndex
                        % fprintf(chainFileUnit                   , self.ProcessID    (i) ...
                                                                % , self.DelRejStage  (i) ...
                                                                % , self.MeanAccRate  (i) ...
                                                                % , self.Adaptation   (i) ...
                                                                % , self.BurninLoc    (i) ...
                                                                % , self.Weight       (i) ...
                                                                % , self.LogFunc      (i) ...
                                                                % , self.State (1:ndim,i) ...
                                                                % ) ;
                    % end
                elseif isVerbose
                    for i = compactStartIndex : compactEndIndex
                        for j = 1 : self.Weight(i)
                            fprintf(chainFileUnit   , chainFileFormat                   ...
                                                    , self.ProcessID    (i)             ...
                                                    , self.DelRejStage  (i)             ...
                                                    , self.MeanAccRate  (i)             ...
                                                    , self.Adaptation   (i)             ...
                                                    , self.BurninLoc    (i)             ...
                                                    , 1                                 ...
                                                    , self.LogFunc      (i)             ...
                                                    , self.State (1:ndim,i)             ...
                                                    ) ;
                        end
                    end
                end
            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end