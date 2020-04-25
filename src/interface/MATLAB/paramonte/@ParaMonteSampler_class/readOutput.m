function outputList = readOutput(self,file,delimiter,fileType)

    self.Err.marginTop = 0;

    if isempty(file)
        if isempty(self.spec.outputFileName)
            self.Err.msg    = "The " + self.methodName + " input simulation specification " + self.objectName + ".spec.outputFileName is not set." ...
                            + "This information is essential, otherwise how could the output files be found?" ...
                            + "All that is needed is the common section of the paths to the output simulation files (including the simulation name)" ...
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
    markovChainRequested = false;
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
        prop="outputListName"; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
        self.outputListName = outputList;
        if ~self.mpiEnabled
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;
            self.Err.msg    = "The processed " + fileType + " file(s) are now stored in the newly-created component """ + outputListFullName + """ of the " + self.methodName + " object. " ...
                            + "For example, to access the contents of the first (or the only) " + fileType + " file, try:" + newline + newline ...
                            + "    " + outputListFullName + "{1}.df" + newline + newline ...
                            + "To access the plotting tools, try:" + newline + newline ...
                            + "    " + outputListFullName + "{1}.plot.<PRESS TAB TO SEE THE LIST OF PLOTS>" + newline + newline ...
                            + "For more information and examples on the usage, visit:" + newline + newline ...
                            + "    <a href=""https://www.cdslab.org/paramonte/"">https://www.cdslab.org/paramonte/</a>";
            self.Err.note();
        end
        outputList = [];
    end

end
