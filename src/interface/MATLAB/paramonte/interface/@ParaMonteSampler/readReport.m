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
%   readReport(file,delimiter)
%
%   Return a list of the contents of a set of ParaMonte simulation(s) output
%   report files whose names begin the user-provided prefix, specified,
%   by the input simulation specification pmpd.spec.outputFileName.
%
%       WARNING
%
%           This method is to be only used for post-processing of the output
%           report file(s) of an already finished simulation. It is NOT meant 
%           to be called by all processes in parallel mode, although it is possible.
%
%   Parameters
%   ----------
%
%       file (optional)
%
%           A string representing the path to the report file with the
%           default value of []. The path only needs to uniquely identify
%           the name of the simulation to which the report file belongs.
%           For example, specifying "./mydir/mysim" as input will lead to
%           a search for a file that begins with "mysim" and ends with
%           "_report.txt" inside the directory "./mydir/".
%
%           If there are multiple files with such name, then all of them
%           will be read and returned as a list.
%
%           If this input argument is not provided by the user, the
%           value of the object's `spec` attribute `outputFileName`
%           will be used instead.
%
%           If the specified path is a URL, the file will be downloaded 
%           as a temporary file to the local system and its contents will 
%           be parsed and the file will be subsequently removed.
%
%           If no input is specified via any of the possible routes, 
%           the method will search for any possible candidate file 
%           with the appropriate suffix in the current working directory.
%
%           Example usage:
%
%               pmpd.readReport();
%
%           or,
%
%               pmpd.readReport("./out/test_run_");
%
%           or,
%
%               pmpd.spec.outputFileName = "./out/test_run_";
%               pmpd.readReport();
%
%           The last two of the above examples are equivalent.
%
%   Returns
%   -------
%
%       reportList (optional)
%
%           A cell array of objects, each of which corresponds to the contents
%           of a unique report file. Each object has the following components:
%
%               file
%
%                   The full absolute path to the report file.
%
%               ndim
%
%                   The number of dimensions of the domain of the objective 
%                   function for which the report was generated.
%
%               contents
%
%                   A char vector containing the entire contents of the report file.
%
%               lineList
%
%                   A one-dimensional cell array whose length is the number of 
%                   the report file. Each element of the array contains one line of 
%                   the report file. 
%
%               setup
%
%                   A structure containing information about the simulation setup as
%                   it appears in the report file, as separate components:
%
%                       io
%
%                           Contains information about the input/output of the 
%                           simulation. 
%
%                       library
%
%                           Contains information about the ParaMonte library that 
%                           was used to generate the report file, including the 
%                           ParaMonte Banner, the simulation environment, the 
%                           computing platform, ... .
%
%               spec
%
%                   A string containing the entire simulation specifications that 
%                   appear in the report file. 
%
%               stats
%
%                   A structure containing a variety of post-processing information
%                   about the simulation that has been performed, including 
%                   but not limited to, 
%
%                       -   timing of the simulation. 
%                       -   the number of function calls of the simulation. 
%                       -   extensive information about the parallelization 
%                           of the simulation, if run in parallel. 
%                       -   summary statistical information about the results 
%                           of the simulation.
%
%       NOTE 
%
%           If no output argument is provided, a reportList property will be added
%           to the parent sampler-object to which the method readReport() belongs.
%           return value of the method. Otherwise, the list will be stored in a
%
%   WARNING 
%
%       If the contents of the report is manipulated, there is no 
%       guarantee of the correct functionality of readReport() method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function reportList = readReport(self, file)

    if isempty(self.objectName); self.objectName = inputname(1); end
    callerName = string(mfilename());
    outputType = string(callerName{1}(5:end));
    outputType = string( [ lower(outputType{1}(1)) , outputType{1}(2:end) ] );
    outputName = outputType + "List";
    fileType = "report";

    % read the report file

    if nargin<2
        file = [];
    else
        if ischar(file) || isstring(file)
            file = convertStringsToChars(file);
        else
            self.Err.msg    = callerName + " takes only one input string or char argument (`file`). Correct usage:" + newline + newline ...
                            + "    " + callerName + "(file)" + newline + newline ...
                            + "where `file` is the path, name, or prefix of the file to be read.";
            self.Err.abort();
        end
    end

    [filePathList, iswebfile] = self.getFilePathList(file,outputType);
    filePathListLen = length(filePathList);

    reportListLocal = cell(filePathListLen,1);
    for ifile = filePathListLen:-1:1
        if self.reportEnabled
            self.Err.msg = "processing file: " + filePathList(ifile);
            self.Err.marginTop = 1;
            self.Err.marginBot = 0;
            self.Err.note();
        end
        reportListLocal{ifile} = ReportFileContents ( filePathList(ifile) ...
                                                    , self.methodName ...
                                                    , self.reportEnabled ...
                                                    , self.Err ...
                                                    );
    end

    % output the reportList

    if nargout==0
        outputListFullName = self.objectName + "." + outputName;
        prop=outputName; if ~any(strcmp(properties(self),prop)); self.addprop(prop); end
        self.(outputName) = reportListLocal;
        self.Err.msg    = "The processed " + fileType + " files are now stored in the newly-created component """ + outputListFullName + """ of the " + self.methodName + " object as a cell array. " ...
                        + "For example, to access the contents of the first (or the only) " + fileType + " file, try:";
    elseif nargout==1
        reportList = reportListLocal;
        outputListFullName = "OUTPUT_CELL_ARRAY";
    else
        self.Err.msg    = "The method, " + self.objectName + "." + callerName + "(file)" ...
                        + "optionally outputs one variable (" + outputName + ") or nothing. If the latter is chosen by the user " ...
                        + "(that is, no output is provivded to the method, " + self.objectName + "." + callerName + "), then the output " + outputName + ...
                        + " will be instead added as a component of the " + self.object + " object.";
        self.Err.abort();
    end

    if self.reportEnabled
        self.Err.marginTop = 1;
        self.Err.marginBot = 1;
        self.Err.msg    = self.Err.msg ...
                        + newline + newline ...
                        + "    " + outputListFullName + "{1}.contents" + newline ...
                        + newline ...
                        + "To access the simulation statistics and information, examine the contents of the" + newline ...
                        + "components of the following structures:" + newline ...
                        + newline ...
                        + "    " + outputListFullName + "{1}.setup             % to get information about the simulation setup." + newline ...
                        + "    " + outputListFullName + "{1}.lineList          % to get the list of lines of the report file." + newline ...
                        + "    " + outputListFullName + "{1}.stats.time        % to get the timing information of the simulation." + newline ...
                        + "    " + outputListFullName + "{1}.stats.chain       % to get the statistics of the simulation output sample." + newline ...
                        + "    " + outputListFullName + "{1}.stats.numFuncCall % to get information about the number of function calls." + newline ...
                        + "    " + outputListFullName + "{1}.stats.parallelism % to get information about the simulation parallelism." + newline ...
                        + "    " + outputListFullName + "{1}.spec              % to get the simulation specification in the report file." + newline ...
                        + newline ...
                        + "For more information and examples on the usage, visit:" + newline + newline ...
                        + "    " + href(self.website.home.url);
        self.Err.note();
    end

    if iswebfile
        self.delFile( filePathList(1) ... file path
                    , "the temporarily-downloaded " + fileType + " file" ... desc
                    );
    end

end