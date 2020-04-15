function fileList = getFileList(self,file,fileType)

    suffix = "_" + fileType + ".txt";
    propertyName = self.objectName + ".spec.outputFileName";

    if exist(file,'file') % check if the input path is a full path to a file

        fileList = [file];
        pattern = file;
        if ~endsWith(file,suffix)
            self.Err.msg    = "The name of the input file:" + newline + newline ...
                            + "    " + file + newline + newline ...
                            + "does not end with the expected suffix '" + suffix + "' for a " + fileType + " file type.";
            self.Err.warn();
        end

    elseif exist(file,'dir') % ensure the input path is not a directory

        self.Err.msg    = propertyName + " = """ + file + """ cannot point to a directory." + newline ...
                        + "Provide a string as the value of file that points to a unique " + fileType + " file or" + newline ...
                        + "to the unique name (including path) of a simulation name that is shared among its output files.";
        self.Err.abort();

    else % search for files matching the input pattern

        if endsWith(file,"*")
            pattern = file;
        else
            pattern = file + "*";
        end

        dirList = dir(pattern);

        if isempty(dirList)
            
        end

        counter = 0;
        fileList = [];
        for ifile = 1:length(dirList)
            if contains(dirList(ifile).name,suffix)
                counter = counter + 1;
                fileList = [ fileList , fullfile( string(dirList(ifile).folder) , string(dirList(ifile).name) ) ];
            end
        end

        if isempty(fileList)
            self.Err.msg    = "Failed to detect any " + fileType + " files with the requested pattern:" + newline + newline ...
                            + "    " + pattern + newline + newline ...
                            + "Provide a string, as the value of the simulation specification " + propertyName + ", that either" + newline + newline ...
                            + "    - points to one or more " + fileType + " files, or," + newline ...
                            + "    - represents the unique name of a " + self.methodName + " simulation." + newline ...
                            + "      This unique-name is the common prefix in the names of" + newline ...
                            + "      the output files of a given " + self.methodName + " simulation.";
            self.Err.abort();
        %else
        %    pattern = pattern + suffix;
        end

    end

    if ~self.mpiEnabled
        self.Err.msg = string(length(fileList)) + " files detected matching the pattern: """ +  strrep(pattern,"\","\\") + """";
        self.Err.note();
    end