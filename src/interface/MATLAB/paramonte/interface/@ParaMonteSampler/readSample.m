function [varargout] = readSample(self,varargin)
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
%       component of the ParaDRAM object named 'sampleList'.
%

    if isempty(self.objectName); self.objectName = inputname(1); end
    callerName = string(mfilename());
    chainType = string(callerName{1}(5:end));
    chainType = string( [ lower(chainType{1}(1)) , chainType{1}(2:end) ] );
    output = chainType + "List";

    if nargout==0
        %if ~any(strcmp(properties(self),output)); self.addprop(output); end
        %self.(output) = self.readOutput(callerName,varargin{:});
        self.readOutput(callerName,varargin{:});
    elseif nargout==1
        %eval(output+" = self.readOutput(callerName,varargin{:})");
        varargout{1} = self.readOutput(callerName,varargin{:});
    else
        self.Err.msg    = "The method, " + self.objectName + "." + callerName + "(file,delimiter)" ...
                        + "optionally outputs one variable (" + output + ") or nothing. If the latter is chosen by the user " ...
                        + "(that is, no output is provivded to the method, " + self.objectName + "." + callerName + "), then the output " + output + ...
                        + " will be instead added as a component of the " + self.object + " object.";
        self.Err.abort();
    end

end