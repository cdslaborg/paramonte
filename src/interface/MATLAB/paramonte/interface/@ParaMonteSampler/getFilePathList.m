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

function [filePathList, iswebfile] = getFilePathList(self,file,fileType)

    self.Err.marginTop = 1;
    self.Err.marginBot = 0;

    iswebfile = false;
    suffix = "_" + fileType + ".txt";
    propertyName = self.objectName + ".spec.outputFileName";

    if isempty(file)
        if isempty(self.spec.outputFileName)
            if ~isempty(self.inputFile)
                self.Err.msg    = self.Err.msg + newline ...
                                + "Apparently, you have specified an input file for the simulation via the attribute """ + self.objectName + ".inputFile""." + newline ...
                                + "Extract the value of outputFileName from this file, and assign it to the simulation specification """ + self.objectName + ".spec.outputFileName"".";
                self.Err.abort();
            end
            if self.reportEnabled
                self.Err.msg    = "The " + self.methodName + " input simulation specification " + self.objectName + ".spec.outputFileName is not set. " ...
                                + "This information is essential, otherwise how could the output files be found? " ...
                                + "All that is needed is the common section of the paths to the output simulation files (including the simulation name) " ...
                                + "or simply, the path to the " + fileType + " file. For now, " + self.methodName + " will search the current working " ...
                                + "for potential " + fileType + " files..." ...
                                ;
                self.Err.warn();
            end
            file = "./";
        else
            file = self.spec.outputFileName;
        end
    end

    filePath = string(getFullPath(convertStringsToChars(file),"lean"));

    if isfile(filePath) % check if the input path is a full path to a file

        filePathList = [filePath];
        pattern = filePath;
        if self.reportEnabled && ~endsWith(filePath,suffix)
            self.Err.msg    = "The name of the input file:" + newline + newline ...
                            + "    """ + filePath + """" + newline + newline ...
                            + "does not end with the expected suffix '" + suffix + "' for a " + fileType + " file type.";
            self.Err.warn();
        end

    else % search for files matching the input pattern

        if isfolder(filePath) % check if the input path is a directory
            if ispc; slash = string('\'); else; slash = "/"; end
            if ~endsWith(filePath,slash); filePath = filePath + slash; end
        end

        if endsWith(filePath,"*")
            pattern = filePath;
        else
            pattern = filePath + "*";
        end

        dirList = dir(pattern);

        counter = 0;
        filePathList = [];
        for ifile = 1:length(dirList)
            if contains(dirList(ifile).name,suffix)
                counter = counter + 1;
                filePathModified = fullfile( string(dirList(ifile).folder) , string(dirList(ifile).name) );
                %filePathModified = string( strrep(filePathModified,'\','\\') );
                filePathList = [ filePathList , filePathModified ];
            end
        end

        if isempty(filePathList) % XXX this may be improved in the future to 
            try
                iswebfile = isurl(file); % check if the input path is a url
            catch
                iswebfile = false;
            end
            if iswebfile
                try
                    filePathList = [ string(websave(genRandFileName(),file)) ];
                    pattern = file;
                catch
                    iswebfile = false;
                    warning( newline + "Failed to read data from the URL: " + file + newline );
                end
            end
        end

        if isempty(filePathList)
            self.Err.msg    = "Failed to detect any " + fileType + " files with the requested pattern:" + newline ...
                            + newline ...
                            + "    " + pattern + newline ...
                            + newline ...
                            + "Provide a string as the value of the simulation specification, " + propertyName + ", that either" + newline ...
                            + newline ...
                            + "    1.   points to one or more " + fileType + " files, or," + newline ...
                            + "    2.   represents the unique name of a " + self.methodName + " simulation." + newline ...
                            + "         This unique-name is the common prefix in the names of" + newline ...
                            + "         the output files of a given " + self.methodName + " simulation." + newline ...
                            + newline ...
                            + "         Note that binary files ending with the suffix ""*.bin"" cannot be parsed." ...
                            ;
            if contains(fileType,"restart")
            self.Err.msg    = self.Err.msg + newline ...
                            + "         You can set the restartFileFormat simulation specification of the sampler to " + newline ...
                            + newline ...
                            + "             " + self.objectName + ".spec.restartFileFormat = ""ascii"";" + newline ...
                            + newline ...
                            + "         to generate ASCII-format MATLAB-readable " + fileType + " files." ...
                            ;
            elseif contains(fileType,"chain")
            self.Err.msg    = self.Err.msg + newline ...
                            + "         You can set the chainFileFormat simulation specification of the sampler to " + newline ...
                            + newline ...
                            + "             " + self.objectName + ".spec.chainFileFormat = ""ascii"";" + newline ...
                            + newline ...
                            + "         to generate ASCII-format MATLAB-readable " + fileType + " files." ...
                            ;
            end
            self.Err.abort();
        %else
        %    pattern = pattern + suffix;
        end

    end

    if self.reportEnabled
        self.Err.msg = string(length(filePathList)) + " files detected matching the pattern: """ +  pattern + """";
        self.Err.note();
    end