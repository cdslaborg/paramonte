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
%
%   This is an internal ParaMonteSampler method and 
%   it is not supposed to be accessible to the users.
%
function varargout = readTabular(self, callerName, file, delimiter) % callerName, file, delimiter

    if isempty(self.objectName); self.objectName = inputname(1); end

    self.Err.marginTop = 0;

    markovChainRequested = false;
    if strcmp(callerName,"readChain")
        fileType = "chain";
    elseif strcmp(callerName,"readSample")
        fileType = "sample";
    elseif strcmp(callerName,"readMarkovChain")
        fileType = "chain";
        markovChainRequested = true;
    elseif strcmp(callerName,"readProgress")
        fileType = "progress";
    end

    if isempty(delimiter)
        if isempty(self.spec.outputDelimiter)
            delimiter = ",";
            if self.reportEnabled
                self.Err.msg    = "The " + self.methodName + " input simulation specification `" + self.objectName + ".spec.outputDelimiter` is not set." + newline ...
                                + "This information is essential for successful reading of the requested " + fileType + " file(s)." + newline ...
                                + "Proceeding with the default assumption of comma-delimited " + fileType + " file contents...";
                self.Err.marginTop = 1;
                self.Err.marginBot = 0;
                self.Err.warn();
            end
        else
            delimiter = string(self.spec.outputDelimiter);
        end
    else
        if ischar(delimiter) || isstring(delimiter)
            delimiter = string(delimiter);
        else
            self.Err.msg    = callerName + " takes only two input string arguments (file, delimiter). Correct usage:" + newline + newline ...
                            + "    " + callerName + "(file, delimiter)" + newline + newline ...
                            + "where `file` is the name of the file to be read, `delimiter` is the delimiter used in the file.";
            self.Err.marginTop = 1;
            self.Err.marginBot = 1;
            self.Err.abort();
        end
    end

    [filePathList, iswebfile] = self.getFilePathList(file,fileType);
    filePathListLen = length(filePathList);

    outputList = cell(filePathListLen,1);
    for ifile = filePathListLen:-1:1
        if self.reportEnabled
            self.Err.msg = "processing file: """ + filePathList(ifile) + """";
            self.Err.marginTop = 1;
            self.Err.marginBot = 0;
            self.Err.note();
        end
        outputList{ifile} = TabularFileContents ( filePathList(ifile) ...
                                                , fileType ...
                                                , delimiter ...
                                                , self.methodName ...
                                                , self.reportEnabled ...
                                                , markovChainRequested ...
                                                , self.Err ...
                                                );
    end

    if nargout==0
        dummy = fileType; if contains(lower(callerName),"markov"); dummy = "markovChain"; end
        outputListName = dummy + "List"; 
        outputListFullName = self.objectName + "." + outputListName;
        prop=outputListName; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
        self.(outputListName) = outputList;
        self.Err.msg    = "The processed " + fileType + " files are now stored in the newly-created component """ + outputListFullName + """ of the " + self.methodName + " object as a cell array. " ...
                        + "For example, to access the contents of the first (or the only) " + fileType + " file, try:";
    elseif nargout==1
        varargout{1} = outputList;
        outputListFullName = "OUTPUT_CELL_ARRAY";
        self.Err.msg    = "The processed " + fileType + " files are now stored in the output variable as a cell array. " ...
                        + "For example, to access the contents of the first (or the only) " + fileType + " file stored in an output variable named " ...
                        + outputListFullName + ", try:";
    end

    if self.reportEnabled
        if strcmp(fileType,"progress")
            specials = "";
        else
            specials    = "    " + outputListFullName + "{1}.plot.line3.make()        % to make 3D line plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.scatter3.make()     % to make 3D scatter plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.lineScatter3.make() % to make 3D line-scatter plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.contour3.make()     % to make 3D kernel-density contour plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.contourf.make()     % to make 2D kernel-density filled-contour plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.contour.make()      % to make 2D kernel-density plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.histogram2.make()   % to make 2D histograms." + newline ...
                        + "    " + outputListFullName + "{1}.plot.histogram.make()    % to make 1D histograms." + newline ...
                        + "    " + outputListFullName + "{1}.plot.grid.make()         % to make GridPlot" + newline ...
                        ;
        end
        self.Err.marginTop = 1;
        self.Err.marginBot = 1;
        self.Err.msg    = self.Err.msg + newline ...
                        + newline ...
                        + "    " + outputListFullName + "{1}.df" + newline ...
                        + newline ...
                        + "To access the plotting tools, try:" + newline ...
                        + newline ...
                        + "    " + outputListFullName + "{1}.plot.<PRESS TAB TO SEE THE LIST OF PLOTS>" + newline ...
                        + newline ...
                        + "For example," + newline ...
                        + newline ...
                        + "    " + outputListFullName + "{1}.plot.line.make()         % to make 2D line plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.scatter.make()      % to make 2D scatter plots." + newline ...
                        + "    " + outputListFullName + "{1}.plot.lineScatter.make()  % to make 2D line-scatter plots." + newline ...
                        + specials ...
                        + newline ...
                        + "To plot or inspect the variable autocorrelations or the correlation/covariance matrices, try:" + newline ...
                        + newline ...
                        + "    " + outputListFullName + "{1}.stats.<PRESS TAB TO SEE THE LIST OF COMPONENTS>" + newline ...
                        + newline ...
                        + "For more information and examples on the usage, visit:" + newline ...
                        + newline ...
                        + "    " + href(self.website.home.url);
        self.Err.note();
    end

    if iswebfile
        self.delFile( filePathList(1) ... file path
                    , "the temporarily-downloaded " + fileType + " file" ... desc
                    );
    end

end
