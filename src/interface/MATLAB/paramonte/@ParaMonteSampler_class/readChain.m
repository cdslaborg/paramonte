function chainList = readChain(self,varargin)
%   Return a list of the contents of a set of ParaDRAM output
%   chain files whose names begin the user-provided prefix, specified,
%   by the input simulation specification pmpd.spec.outputFileName
%   This method is to be only used for postprocessing of the output
%   chain file(s) of an already finished ParaDRAM simulation.
%   It is not meant to be called by all processes in parallel mode,
%   although it is possible.
%
%   Parameters
%   ----------
%       file
%           A string representing the path to the chain file with
%           the default value of None.
%           The path only needs to uniquely identify the simulation
%           to which the chain file belongs. For example, specifying
%           "./mydir/mysim" as input will lead to a search for a file
%           that begins with "mysim" and ends with "_chain.txt"
%           inside the directory "./mydir/". If there are multiple
%           files with such name, then all of them will be read
%           and returned as a list.
%           If this input argument is not provided by the user, the
%           value of the object attribute outputFileName
%           will be used instead. At least one of the two mentioned
%           routes must provide the path to the chain file otherwise,
%           this method will break by calling sys.exit().
%       delimiter
%           Optional input string representing the delimiter used in the
%           output chain file. If it is not provided as input argument,
%           the value of the corresponding object attribute outputDelimiter
%           will be used instead. If none of the two are available,
%           the default comma delimiter "," will be assumed and used.
%       parseContents
%           If set to True, the contents of the file will be parsed and
%           stored in a component of the object named 'contents'.
%           The default value is True.
%       renabled
%           If set to False, the contents of the file(s) will be stored as a
%           list in a (new) component of the ParaDRAM object named 'chainList'
%           and None will be the return value of the method.
%           If set to True, the reverse will done.
%           The default value is False.
%
%   Returns
%   -------
%       a list of objects, each of which has the following properties:
%           file
%               full absolute path to the chain file.
%           delimiter
%               the delimiter used in the chain file.
%           ndim
%               number of dimensions of the domain of the objective function
%               from which the chain has been drawn.
%           count
%               the number of unique (weighted) points in the chain file.
%               This is essentially the number of rows in the chain file
%               minus one (representing the header line).
%           df
%               the contents of the chain file in the form of
%               a pandas-library DataFrame (hence called 'df').
%           dynamic attributes:
%               corresponding to each column in the chain file, a property
%               with the same name as the column header is also created
%               for the object which contains the data stored in that column
%               of the chain file.
%       If renabled = True, the list of objects will be returned as the
%       return value of the method. Otherwise, the list will be stored in a
%       component of the ParaDRAM object named 'chainList'.
%

    if isempty(self.objectName)
        self.objectName = inputname(1);
    end

    %p = inputParser;
    %renabledValueDefault = false;
    %renabledValueValid = {true,false};
    %renabledCheck = @(x) any(validatestring(x,renabledValueValid));
    %addOptional(p,'renabled',renabledValueDefault,renabledCheck);
    %fileValueDefault = false;
    %fileValueValid = {char,string};
    %fileCheck = @(x) any(validatestring(x,fileValueValid));
    %addOptional(p,'file',fileValueDefault,fileCheck);
    %p.KeepUnmatched = true;
    %parse(p,file,varargin{:})
    %parse(p,renabled,varargin{:})

    file = [];
    delimiter = [];
    renabled = [];
    scriptName = mfilename();
    nargin
    if nargin>1
        file = varargin{1};
    end
    if nargin>2
        delimiter = varargin{2};
    end
    if nargin>3
        renabled = varargin{3};
    end
    if nargin>4
        self.Err.msg    = "More than 3 input arguments detected. Correct usage:" + newline + newline ...
                        + "    " + self.objectName + "." + scriptName + "(file,delimiter,renabled)" + newline + newline ...
                        + "where ""file"" is the name of the file to be read, delimiter is the delimiter used in the file," + newline ...
                        + "and renabed, if set to true, will return the file contents as a component of the " + self.methodName + " " + self.object + " instance.";
        self.Err.abort();
    end

    if isempty(renabled)
        self.readOutput(file,delimiter,"chain");
    else
        outputList = self.readOutput(file,delimiter,"chain");
    end

end