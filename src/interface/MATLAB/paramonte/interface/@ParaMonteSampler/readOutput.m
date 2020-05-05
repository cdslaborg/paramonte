function varargout = readOutput(self,varargin) % callerName,file,delimiter

    if isempty(self.objectName); self.objectName = inputname(1); end

    self.Err.marginTop = 0;

    file = [];
    delimiter = [];
    errorOccurred = false;
    markovChainRequested = false;
    if nargin==1; errorOccurred = true; end
    if nargin>1
        callerName = varargin{1};
        if strcmp(callerName,"readChain")
            fileType = "chain";
        elseif strcmp(callerName,"readSample")
            fileType = "sample";
        elseif strcmp(callerName,"readMarkovChain")
            fileType = "chain";
            markovChainRequested = true;
        end
    end
    if nargin>2; file = varargin{2}; end
    if nargin>3; delimiter = varargin{3}; end
    if nargin>4 || errorOccurred
        self.Err.msg    = callerName + "takes only two input arguments (file, delimiter). Correct usage:" + newline + newline ...
                        + "    " + callerName + "(file,delimiter)" + newline + newline ...
                        + "where ""file"" is the name of the file to be read, delimiter is the delimiter used in the file,";
        self.Err.abort();
    end

    if isempty(file)
        if isempty(self.spec.outputFileName)
            self.Err.msg    = "The " + self.methodName + " input simulation specification " + self.objectName + ".spec.outputFileName is not set. " ...
                            + "This information is essential, otherwise how could the output files be found? " ...
                            + "All that is needed is the common section of the paths to the output simulation files (including the simulation name) " ...
                            + "or simply, the path to the " + fileType + " file.";
            if ~isempty(self.inputFile)
                self.Err.msg    = self.Err.msg + newline ...
                                + "Apparently, you have specified an input file for the simulation via the attribute """ + self.objectName + ".inputFile""." + newline ...
                                + "Extract the value of outputFileName from this file, and assign it to the simulation specification """ + self.objectName + ".spec.outputFileName"".";
            end
            self.Err.abort();
        else
            file = string(self.spec.outputFileName);
        end
    else
        file = string(file);
    end

    if isempty(delimiter)
        if isempty(self.spec.outputDelimiter)
            delimiter = ',';
            if ~self.mpiEnabled
                self.Err.msg    = "The " + self.methodName + " input simulation specification " + self.objectName + ".spec.outputDelimiter is not set." + newline ...
                                + "This information is essential for successful reading of the requested " + fileType + " file(s)." + newline ...
                                + "Proceeding with the default assumption of comma-delimited " + fileType + " file contents...";
                self.Err.warn();
            end
        else
            delimiter = string(self.spec.outputDelimiter);
        end
    else
        delimiter = string(delimiter);
    end

    fileList = self.getFileList(file,fileType);
    fileListLen = length(fileList);

    outputList = cell(fileListLen,1);
    for ifile = fileListLen:-1:1
        filePathModified = string( strrep(fileList(ifile),'\','\\') );
        if ~self.mpiEnabled
            self.Err.msg = "processing file: " + filePathModified;
            self.Err.marginTop = 1;
            self.Err.marginBot = 0;
            self.Err.note();
        end
        outputList{ifile} = OutputFileContents  ( filePathModified ...
                                                , fileType ...
                                                , delimiter ...
                                                , self.methodName ...
                                                , self.mpiEnabled ...
                                                , markovChainRequested ...
                                                , self.Err ...
                                                );
    end

    if nargout==0
        outputListName = fileType + "List";
        outputListFullName = self.objectName + "." + outputListName;
        prop=outputListName; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
        self.(outputListName) = outputList;
        self.Err.msg    = "The processed " + fileType + " file(s) are now stored in the newly-created component """ + outputListFullName + """ of the " + self.methodName + " object as a cell array. " ...
                        + "For example, to access the contents of the first (or the only) " + fileType + " file, try:";
    elseif nargout==1
        varargout{1} = outputList;
        outputListFullName = "OUTPUT_CELL_ARRAY";
        self.Err.msg    = "The processed " + fileType + " file(s) are now stored in the output variable as a cell array. " ...
                        + "For example, to access the contents of the first (or the only) " + fileType + " file stored in an output variable named " ...
                        + outputListFullName + ", try:";
    end

    if ~self.mpiEnabled
        self.Err.marginTop = 1;
        self.Err.marginBot = 1;
        self.Err.msg    = self.Err.msg ...
                        + newline + newline ...
                        + "    " + outputListFullName + "{1}.df" + newline + newline ...
                        + "To access the plotting tools, try:" + newline + newline ...
                        + "    " + outputListFullName + "{1}.plot.<PRESS TAB TO SEE THE LIST OF PLOTS>" + newline + newline ...
                        + "To plot or inspect the variable autocorrelations or the correlation/covariance matrices, try:" + newline + newline ...
                        + "    " + outputListFullName + "{1}.stats.<PRESS TAB TO SEE THE LIST OF COMPONENTS>" + newline + newline ...
                        + "For more information and examples on the usage, visit:" + newline + newline ...
                        + "    <a href=""https://www.cdslab.org/paramonte/"">https://www.cdslab.org/paramonte/</a>";
        self.Err.note();
    end

end
