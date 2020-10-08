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

classdef ParaMonte_class < dynamicprops % handle % dynamicprops is needed for adding postprocessing properties 

    properties (Constant, Access = protected, Hidden)
        CLASS_NAME = "@ParaMonte_class"
    end

    properties (Access = public, Hidden)
        nd                      = []
        name                    = []
        brand                   = []
        date                    = []
        version                 = []
        isDryRun                = []
        isFreshRun              = []
        OS                      = []
        Err                     = []
        Image                   = []
        SpecBase                = []
        SystemInfo              = []
        Timer                   = []
        LogFile                 = []
        TimeFile                = []
        ChainFile               = []
        SampleFile              = []
        RestartFile             = []
        Decor                   = []
        website                 = []
        platform                = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public, Hidden)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = ParaMonte_class(platform,website)

            self.nd                         = [];
            self.name                       = [];
            self.brand                      = [];
            self.date                       = [];
            self.version                    = [];
            self.isDryRun                   = [];
            self.isFreshRun                 = [];
            self.OS                         = OS_class();
            self.Err                        = Err_class();
            self.Image                      = struct( 'id', []          ...
                                                    , 'count', []       ...
                                                    , 'isFirst', []     ...
                                                    , 'isNotFirst', []  ...
                                                    , 'isMaster', []    ...
                                                    , 'isNotMaster', [] ...
                                                    , 'name', []        ...
                                                    ) ;


            self.Timer                      = Timer_class();
            self.LogFile                    = struct('unit', [], 'exists', false, 'suffix', "report");
            self.TimeFile                   = struct('unit', [], 'exists', false, 'suffix', "progress");
            self.ChainFile                  = struct('unit', [], 'exists', false, 'suffix', "chain");
            self.SampleFile                 = struct('unit', [], 'exists', false, 'suffix', "sample");
            self.RestartFile                = struct('unit', [], 'exists', false, 'suffix', "restart", 'counter', 0);
            self.Decor                      = Decoration_class([],[],[],[]);
            self.website = website;
            self.platform = platform;

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function setupParaMonte(self, nd, name, date, version)

            FUNCTION_NAME           = self.CLASS_NAME + "@setupParaMonte()";

            self.nd.val             = nd;

            self.LogFile.unit       = 1;            % temporarily set the report file to stdout.
            self.Decor              = Decoration_class([],[],[],[]);  % initialize the TAB character and decoration symbol to the default values.

            self.name               = name;
            self.date               = date;
            self.brand              = Decoration_class.TAB + Decoration_class.TAB + self.name;
            self.version            = "Version " + version;

            self.Image.id           = 1;
            self.Image.count        = 1;

            self.Image.name         = "@process(" + num2str(self.Image.id) + ")";
            self.Image.isFirst      = self.Image.id == 1;
            self.Image.isNotFirst   = ~self.Image.isFirst;
            self.Image.isMaster     = true;    % ATTN in Parallel mode: set to false initially and change later on, depending on the requested type of parallelism
            self.Image.isNotMaster  = ~self.Image.isMaster;

            % setup formatting variables

            self.nd.str             = num2str(self.nd.val);

            % determine OS

            self.OS.queryOS();
            if self.OS.Err.occurred
                self.Err            = self.OS.Err;
                self.Err.msg        = FUNCTION_NAME + ": Error occurred while querying OS type." + Constants.NLC + self.Err.msg;
                self.Err.prefix     = self.brand;
                self.Err.newLine    = newline;
                self.Err.outputUnit = self.LogFile.unit;
                self.abort();
            end

            % get system info by all images

            self.SystemInfo = SystemInfo_class(self.OS.isWindows, self.platform.systemInfoFilePath);

            if self.SystemInfo.Err.occurred
                self.Err            = self.SystemInfo.Err;
                self.Err.msg        = FUNCTION_NAME + ": Error occurred while collecting system info." + Constants.NLC + self.Err.msg;
                self.Err.prefix     = self.brand;
                self.Err.newLine    = newline;
                self.Err.outputUnit = self.LogFile.unit;
                self.abort();
            end

            if self.Image.isFirst    % blockSplashByFirstImage
                self.addSplashScreen();
                self.noteUserAboutEnvSetup();
            end

            function index_val = index(str, sub_str)
                index_val   = strfind(str, sub_str);
                if isempty(index_val), index_val = 0; else, index_val = index_val(1); end
            end

            % Set the default and null values for ParaMonte SpecBase

            self.SpecBase = SpecBase_class(self.nd.val, self.name);

        end % function setupParaMonte

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function addSplashScreen(self)
            self.Decor.text = " " + newline + newline...
                            ...%self%name // "\n" // &
                            + "ParaMonte" + newline...
                            + self.version + newline...
                            + self.date + newline...
                            + newline...
                            + "Department of Physics" + newline...
                            + "Computational & Data Science Lab" + newline...
                            + "Data Science Program, College of Science" + newline...
                            + "The University of Texas at Arlington" + newline...
                            + newline...
                            + "originally developed at" + newline...
                            + newline...
                            + "Multiscale Modeling Group" + newline...
                            + "Center for Computational Oncology (CCO)" + newline...
                            + "Oden Institute for Computational Engineering and Sciences" + newline...
                            + "Department of Aerospace Engineering and Engineering Mechanics" + newline...
                            + "Department of Neurology, Dell-Seton Medical School" + newline...
                            + "Department of Biomedical Engineering" + newline...
                            + "The University of Texas at Austin" + newline...
                            + newline...
                            + "For questions and further information, please contact:" + newline...
                            + newline...
                            + "Shashank D. Kumbhare" + ", "...
                            + "Amir Shahmoradi" + newline...
                            + newline...
                            + "shashank.kumbhare@mavs.uta.edu" + newline...
                            + "shashank.kumbhare8@gmail.com" + newline...
                            + "shahmoradi@utexas.edu" + newline...
                            + newline...
                            + "cdslab.org/pm" + newline...
                            + newline...
                            + "https://www.cdslab.org/paramonte/" + newline...
                            + newline...
                            ;

            self.Decor.writeDecoratedText(self.Decor.text   ...
                                        , "*"               ...
                                        , 132               ...
                                        , 4                 ...
                                        , 2                 ...
                                        , 1                 ...
                                        , 2                 ...
                                        , self.LogFile.unit ...
                                        , "\n"              ...
                                        );

        end % function addSplashScreen

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function addCompilerPlatformInfo(self)

            % report the ParaMonte interpretor version and options

           %self.Decor.writeDecoratedText   (  " " + newline + self.name + " library compiler version" + newline ...
            self.Decor.writeDecoratedText   (  " " + newline + " ParaMonte library interface specifications" + newline ...
                                            , "*"                                               ...
                                            , 132                                               ...
                                            , 4                                                 ...
                                            , 1                                                 ...
                                            , 2                                                 ...
                                            , 1                                                 ...
                                            , self.LogFile.unit                                 ...
                                            , "\n"                                              ...
                                            );

            fprintf(self.LogFile.unit, "MATLAB - " + version + newline);

            self.Decor.writeDecoratedText   ( " " + newline + "Runtime platform specifications" + newline   ...
                                            , "*"                                       ...
                                            , 132                                       ...
                                            , 4                                         ...
                                            , 1                                         ...
                                            , 2                                         ...
                                            , 1                                         ...
                                            , self.LogFile.unit                         ...
                                            , "\n"                                      ...
                                            );

            fprintf(self.LogFile.unit, self.SystemInfo.info);

            self.Decor.write(self.LogFile.unit, [], [], [], []);

        end % function addCompilerPlatformInfo

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function noteUserAboutEnvSetup(self)
            self.Decor.writeDecoratedText   ( " " + newline + "Setting up " + self.name + " simulation environment" + newline   ...
                                            , []                                                                                ...
                                            , []                                                                                ...
                                            , []                                                                                ...
                                            , []                                                                                ...
                                            , 1                                                                                 ...
                                            , 1                                                                                 ...
                                            , self.LogFile.unit                                                                 ...
                                            , newline                                                                           ...
                                            );
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function warnUserAboutMissingNamelist(self, prefix, name, namelist, outputUnit)
            msg =   "No namelist group of variables named " + namelist + " was detected in user's input file for " + name + " options." + newline   ...
                    + "All " + name + " options will be assigned appropriate default values."...
                    ;
            self.Err.reset();
            self.Err.prefix      = prefix;
            self.Err.outputUnit  = outputUnit;
            self.Err.msg         = msg;
            self.Err.warn();
            if outputUnit ~= 1
                Err.outputUnit  = 1;
                Err.warn();
            end
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function setupOutputFiles(self)

            FUNCTION_NAME = self.CLASS_NAME + "@setupOutputFiles()";

            if self.SpecBase.outputFileName.original == self.SpecBase.outputFileName.def
                msg = "No user-input filename prefix for " + self.name + " output files detected." + Constants.NLC            ...
                    + "Generating appropriate filenames for " + self.name + " output files from the current date and time..." ...
                    ;
            else
                msg = "Variable outputFileName detected among the input variables to " + self.name + ":" + Constants.NLC + self.SpecBase.outputFileName.original;
            end

            self.SpecBase.outputFileName.queryPath(self.SpecBase.outputFileName.original, self.OS);
            
            self.Err.prefix         = self.brand;
            self.Err.newLine        = newline;
            self.Err.outputUnit     = self.LogFile.unit;
            self.Err.resetEnabled   = false;

            if self.SpecBase.outputFileName.Err.occurred
                self.Err            = self.SpecBase.outputFileName.Err;
                self.Err.msg        = FUNCTION_NAME + ": Error occurred while attempting to construct outputFileName path type." + Constants.NLC + self.Err.msg;
                self.Err.abort();
            end

            self.SpecBase.outputFileName.namePrefix = self.SpecBase.outputFileName.name + self.SpecBase.outputFileName.ext;

            % get the current working directory

            currentWorkingDir       = pwd;
            currentWorkingDir       = strtrim(currentWorkingDir);

            if isempty(currentWorkingDir)
                self.Err.msg        = FUNCTION_NAME + ":Error occurred while fetching the current working directory via pwd." + Constants.NLC;
                self.Err.abort();
            end

            msg = msg + Constants.NLC + Constants.NLC + "Absolute path to the current working directory:" + Constants.NLC + strrep(currentWorkingDir, '\', '\\');

            if isempty(strtrim(char(self.SpecBase.outputFileName.dir)))
                self.SpecBase.outputFileName.dir = convertCharsToStrings(currentWorkingDir); % + self.SpecBase.outputFileName.slashOS;
                msg = msg + Constants.NLC + Constants.NLC + "All output files will be written to the current working directory:" + Constants.NLC + strrep(self.SpecBase.outputFileName.dir, '\', '\\');
            else
                msg = msg + Constants.NLC + Constants.NLC + "Generating the requested directory for ParaDRAM output files:" + Constants.NLC + strrep(self.SpecBase.outputFileName.dir, '\', '\\');
            end

            % Generate the output files directory:

            if self.Image.isFirst && self.SpecBase.outputFileName.dir ~= pwd
                if ~exist(self.SpecBase.outputFileName.dir, 'dir')
                    [status, self.Err.msg]  = mkdir(self.SpecBase.outputFileName.dir);
                    self.Err.occurred       = ~status;
                    if self.Err.occurred
                        self.Err.msg        = FUNCTION_NAME + ": Error occurred while making directory = '" + self.SpecBase.outputFileName.dir + "'." + Constants.NLC + self.Err.msg;
                        self.Err.abort();
                    end
                end
            end

            if isempty(self.spec.outputFileName)
                self.SpecBase.outputFileName.dir = convertCharsToStrings(currentWorkingDir) + self.SpecBase.outputFileName.slashOS;
            end

            if isempty(strtrim(char(self.SpecBase.outputFileName.namePrefix)))
                msg =  msg + Constants.NLC + Constants.NLC + "No user-input filename prefix for " + self.name + " output files detected." + Constants.NLC...
                    + "Generating appropriate filenames for " + self.name + " output files from the current date and time..."                         ...
                    ;
                self.SpecBase.outputFileName.namePrefix = self.SpecBase.outputFileName.def;
            end

            self.SpecBase.outputFileName.pathPrefix = self.SpecBase.outputFileName.dir + self.SpecBase.outputFileName.slashOS + self.SpecBase.outputFileName.namePrefix;


            % Variable msg will be used down this subroutine, so it should not be changed beyond this point

            msg  =  msg + Constants.NLC + Constants.NLC + self.name + " output files will be prefixed with:" + Constants.NLC + strrep(self.SpecBase.outputFileName.pathPrefix,'\',"\\");

            if self.Image.isFirst
                self.Err.msg        = msg;
                self.Err.note();
            end


            % Generate the output filenames, search for pre-existing runs, and open the report file:

            if self.Image.isFirst
                self.Err.msg        = "Searching for previous runs of " + self.name + "...";
                self.Err.note();
            end

            % this block could be all executed by only the master images

            self.LogFile.Path.ext           = Constants.FILE_EXT.ascii;
            self.TimeFile.Path.ext          = Constants.FILE_EXT.ascii;
            self.ChainFile.Path.ext         = Constants.FILE_EXT.ascii;
            self.SampleFile.Path.ext        = Constants.FILE_EXT.ascii;
            self.RestartFile.Path.ext       = Constants.FILE_EXT.ascii;

            self.RestartFile.Form.value     = "formatted";
            if self.SpecBase.restartFileFormat.isBinary
                self.RestartFile.Form.value = "unformatted";
                self.RestartFile.Path.ext   = Constants.FILE_EXT.binary;
            end

            self.ChainFile.Form.value       = "formatted";
            if self.SpecBase.chainFileFormat.isBinary
                self.ChainFile.Form.value   = "unformatted";
                self.ChainFile.Path.ext     = Constants.FILE_EXT.binary;
            end

            imageID = 1;

            fullOutputFileName              = self.SpecBase.outputFileName.pathPrefix + "_process_" + num2str(imageID) + "_";

            self.LogFile.Path               = Path_class(fullOutputFileName + char(self.LogFile.suffix)        + char(self.LogFile.Path.ext)    , self.OS );
            self.TimeFile.Path              = Path_class(fullOutputFileName + char(self.TimeFile.suffix)       + char(self.TimeFile.Path.ext)   , self.OS );
            self.ChainFile.Path             = Path_class(fullOutputFileName + char(self.ChainFile.suffix)      + char(self.ChainFile.Path.ext)  , self.OS );
            self.SampleFile.Path            = Path_class(fullOutputFileName + char(self.SampleFile.suffix)     + char(self.SampleFile.Path.ext) , self.OS );
            self.RestartFile.Path           = Path_class(fullOutputFileName + char(self.RestartFile.suffix)    + char(self.RestartFile.Path.ext), self.OS );

            self.LogFile.exists     = isfile(char(self.LogFile.Path.original));
%            if self.Err.occurred
%                self.Err.msg        = FUNCTION_NAME + ": Error occurred while inquiring the existence of file='"...
%                                    + strrep(self.LogFile.Path.original, '\', '\\') + "'";
%                self.Err.abort();
%            end
%
            self.SampleFile.exists  = isfile(char(self.SampleFile.Path.original));
%            if self.Err.occurred
%                self.Err.msg        = FUNCTION_NAME + ": Error occurred while inquiring the existence of file='"...
%                                    + strrep(self.SampleFile.Path.original, '\', '\\') + "'";
%                self.Err.abort();
%            end
%
            self.TimeFile.exists    = isfile(char(self.TimeFile.Path.original));
%            if self.Err.occurred
%                self.Err.msg        = FUNCTION_NAME + ": Error occurred while inquiring the existence of file='"...
%                                    + strrep(self.TimeFile.Path.original, '\', '\\') + "'";
%                self.Err.abort();
%            end
%
            self.ChainFile.exists   = isfile(char(self.ChainFile.Path.original));
%            if self.Err.occurred
%                self.Err.msg        = FUNCTION_NAME + ": Error occurred while inquiring the existence of file='"...
%                                    + strrep(self.ChainFile.Path.original, '\', '\\') + "'";
%                self.Err.abort();
%            end
%
            self.RestartFile.exists = isfile(char(self.RestartFile.Path.original));
%            if self.Err.occurred
%                self.Err.msg        = FUNCTION_NAME + ": Error occurred while inquiring the existence of file='"...
%                                    + strrep(self.RestartFile.Path.original, '\', '\\') + "'";
%                self.Err.abort();
%            end

            self.isDryRun = ~self.SpecBase.overwriteRequested.val & ( ... not fresh, if any file exists
                            self.LogFile.exists || self.TimeFile.exists || self.RestartFile.exists || self.ChainFile.exists || self.SampleFile.exists);
            self.isFreshRun = ~self.isDryRun;
%            if self.isFreshRun, self.SpecMCMC.chainSize.val = self.SpecMCMC.chainSize.val + 1; end

            if self.isFreshRun
                if self.Image.isFirst
                    self.Err.msg        = "No pre-existing " + self.name + " run detected." + newline + "Starting a fresh ParaDRAM run...";
                    self.Err.note();
                end
            else
                if self.Image.isFirst
                    self.Err.msg = "Previous run of " + self.name + " detected." + Constants.NLC + "Searching for restart files...";
                    self.Err.note();
                end
                if self.SampleFile.exists % sampling is already complete
                    self.Err.occurred   = true;
                    self.Err.msg        = strrep(FUNCTION_NAME + ": Error occurred. Output sample file detected: " + self.SampleFile.Path.original  ...
                                        + Constants.NLC + self.name + " cannot overwrite an already-completed simulation."                          ...
                                        + Constants.NLC + "Please provide an alternative file name for the new simulation outputs.", '\', '\\')     ...
                                        ;
                elseif self.LogFile.exists && self.TimeFile.exists && self.RestartFile.exists && self.ChainFile.exists % restart mode
                    self.Err.occurred   = false;
                else
                    self.Err.occurred   = true;
                    self.Err.msg = FUNCTION_NAME + ": Error occurred. For a successful simulation restart, all output files are necessary." + Constants.NLC + "List of missing simulation output files:";
                    if ~self.LogFile.exists,        self.Err.msg = strrep(self.Err.msg + Constants.NLC + self.LogFile.Path.original,    '\', '\\');       end
                    if ~self.TimeFile.exists,       self.Err.msg = strrep(self.Err.msg + Constants.NLC + self.TimeFile.Path.original,   '\', '\\');      end
                    if ~self.ChainFile.exists,      self.Err.msg = strrep(self.Err.msg + Constants.NLC + self.ChainFile.Path.original,  '\', '\\');     end
                    if ~self.RestartFile.exists,    self.Err.msg = strrep(self.Err.msg + Constants.NLC + self.RestartFile.Path.original,'\', '\\');   end
                end
                if self.Err.occurred, self.Err.abort(); end
            end


            % open/append the output files:

            if self.Image.isMaster
                if self.isFreshRun
                    workingOn                           = "Generating the output ";
                    self.LogFile.       status          = "w";
                    self.TimeFile.      status          = "w";
                    self.ChainFile.     status          = "w";
                    self.SampleFile.    status          = "w";
                    self.RestartFile.   status          = "w";
                    %self.LogFile.       Position.value  = "asis";
                    %self.TimeFile.      Position.value  = "asis";
                    %self.ChainFile.     Position.value  = "asis";
                    %self.SampleFile.    Position.value  = "asis";
                    %self.RestartFile.   Position.value  = "asis";
                else
                    workingOn                           = "Appending to the existing "  ;
                    self.LogFile.       status          = "r+";
                    self.TimeFile.      status          = "r+";
                    self.ChainFile.     status          = "r+";
                    self.SampleFile.    status          = "w";
                    self.RestartFile.   status          = "r+";
                    self.LogFile.       Position.value  = "a+";
                    %self.TimeFile.      Position.value  = "asis";
                    %self.ChainFile.     Position.value  = "asis";
                    %self.SampleFile.    Position.value  = "asis";
                    %self.RestartFile.   Position.value  = "asis";
                end
            end

            % print the stdout message for generating / appending the output report file

            if self.Image.isFirst     % LogFileListByFirstImage

                % print the stdout message for generating / appending the output report file(s)

                self.Err.msg        = workingOn + char(self.LogFile.suffix) + " file:";
                self.Err.marginBot  = 0;
                self.Err.note();

                % print the the output report file name of the images

                self.Err.msg        = strrep(self.LogFile.Path.original, '\', '\\');
                self.Err.marginTop  = 0;
                self.Err.marginBot  = 0;
                self.Err.note();

                self.Err.msg        = "Please see the output " + self.LogFile.suffix + " and " + self.TimeFile.suffix + " files for further realtime simulation details.";
                self.Err.marginTop  = 3;
                self.Err.marginBot  = 3;
                self.Err.note();

            end

            % open the output files

            self.Err.marginTop      = 1;
            self.Err.marginBot      = 1;
            
            if self.Image.isMaster
                
                % open LogFile
                %if self.LogFile.exists,     opentype = "a+";    else, opentype = "w"; end
                [self.LogFile.unit, self.Err.msg] = fopen(self.LogFile.Path.original, self.LogFile.status);
                self.Err.outputUnit = self.LogFile.unit;
                if self.Err.msg
                    self.Err.msg    = FUNCTION_NAME + ": Error occurred while opening " + self.name + self.LogFile.suffix + " file='" + strrep(self.LogFile.Path.original, '\', '\\') + "'.";
                    self.Err.abort();
                end
                % rewrite the same old stuff to all report files
                if self.isFreshRun
                    self.addSplashScreen();
                    if self.SpecBase.silentModeRequested.isFalse, self.addCompilerPlatformInfo(); end   % this takes about 0.75 seconds to execute on Stampede Login nodes.
                    self.noteUserAboutEnvSetup();
                    self.Err.msg    = msg;
                    self.Err.note();
                end

                % workingOn on TimeFile
                if self.isFreshRun
                    self.Err.msg    = workingOn + self.TimeFile.suffix + " file:" + Constants.NLC + strrep(self.TimeFile.Path.original, '\', '\\');
                    self.Err.note();
                end
                % open TimeFile
                %if self.TimeFile.exists,    opentype = "r+";    else, opentype = "w"; end
                [self.TimeFile.unit, self.Err.msg] = fopen(self.TimeFile.Path.original, self.TimeFile.status);
                if self.Err.msg
                    self.Err.msg    = FUNCTION_NAME + ": Error occurred while opening " + self.name + self.TimeFile.suffix + " file='" + strrep(self.TimeFile.Path.original, '\', '\\') + "'.";
                    self.Err.abort();
                end

                % workingOn on ChainFile
                if self.isFreshRun
                    self.Err.msg    = workingOn + self.ChainFile.suffix + "file:" + Constants.NLC + strrep(self.ChainFile.Path.original, '\', '\\');
                    self.Err.note();
                end
                % open ChainFile
                if self.ChainFile.exists,   opentype = "r+";    else, opentype = "w"; end
                [self.ChainFile.unit, self.Err.msg] = fopen(self.ChainFile.Path.original, self.ChainFile.status);
                if self.Err.msg
                    self.Err.msg    = FUNCTION_NAME + ": Error occurred while opening " + self.name + self.ChainFile.suffix + " file='" + strrep(self.ChainFile.Path.original, '\', '\\') + "'.";
                    self.Err.abort();
                end

                % open RestartFile
                %if self.RestartFile.exists, opentype = "r+";    else, opentype = "w"; end
                [self.RestartFile.unit, self.Err.msg] = fopen(self.RestartFile.Path.original, self.RestartFile.status);
                if self.Err.msg
                    self.Err.msg = FUNCTION_NAME + ": Error occurred while opening " + self.name + self.RestartFile.suffix + " file='" + strrep(self.RestartFile.Path.original, '\', '\\') + "'.";
                    self.Err.abort();
                end

                if self.isFreshRun
                    self.Decor.writeDecoratedText   ( " " + Constants.NLC + self.name + " simulation specifications" + Constants.NLC    ...
                                                    , [], [], [], [], 1, 1                                                              ...
                                                    , self.LogFile.unit                                                                 ...
                                                    , Constants.NLC                                                                     ...
                                                    ) ;
                end

            end

            % These must be defined for all images, because they may be passed as arguments to the kernel subroutines.

            formatInt   = "%" + self.SpecBase.outputColumnWidth.str + "d";                                                  % For Ex. "%4d"
            formatReal  = "%" + self.SpecBase.outputColumnWidth.str + "." + self.SpecBase.outputRealPrecision.str + "f";    % For Ex. "%10.8f"
            formatStr   = "%" + self.SpecBase.outputColumnWidth.str + "s";                                                  % For Ex. "%12s"
            delim       = self.SpecBase.outputDelimiter.val;                                                                % Delimiter

            %-----------------------------------------------------------------------------------------------------------------------
            % ChainFile - headerFormat
            self.ChainFile.headerFormat     = formatStr;
            for i = 1 : length(self.Chain.ColHeader)-1
                self.ChainFile.headerFormat = self.ChainFile.headerFormat + delim + formatStr;
                if i == length(self.Chain.ColHeader)-1, self.ChainFile.headerFormat = self.ChainFile.headerFormat + "\n"; end
            end
            
            % ChainFile - format
            self.ChainFile.format           = formatInt + delim + formatInt + delim + formatReal + delim + formatReal + delim + formatInt + delim + formatInt + delim + formatReal;
            for i = 1 : length(self.SpecBase.variableNameList.Val)
                self.ChainFile.format       = self.ChainFile.format + delim + formatReal;
                if i == length(self.SpecBase.variableNameList.Val), self.ChainFile.format = self.ChainFile.format + "\n"; end
            end
            %-----------------------------------------------------------------------------------------------------------------------
            % TimeFile - headerFormat
            self.TimeFile.headerFormat      = formatStr + delim + formatStr + delim + formatStr + delim + formatStr + delim + formatStr + delim + formatStr + delim + formatStr + "\n";
            % TimeFile - format
            self.TimeFile.format            = formatInt + delim + formatInt + delim + formatReal + delim + formatReal + delim + formatReal + delim + formatReal + delim + formatReal + "\n";
            %-----------------------------------------------------------------------------------------------------------------------
            % SampleFile - headerFormat
            self.SampleFile.headerFormat    = formatStr;
            for i = 1 : length(self.SpecBase.variableNameList.Val)
                self.SampleFile.headerFormat= self.SampleFile.headerFormat + delim + formatStr;
                if i == length(self.SpecBase.variableNameList.Val), self.SampleFile.headerFormat = self.SampleFile.headerFormat + "\n"; end
            end
            % SampleFile - format
            self.SampleFile.format          = formatReal;
            for i = 1 : length(self.SpecBase.variableNameList.Val)
                self.SampleFile.format  = self.SampleFile.format + delim + formatReal;
                if i == length(self.SpecBase.variableNameList.Val), self.SampleFile.format = self.SampleFile.format + "\n"; end
            end

            %-----------------------------------------------------------------------------------------------------------------------
            % RestartFile - format
            precisionRestart    = "\n" + "%" + "." + self.SpecBase.outputRealPrecision.str + "f";
            precisionRestart    = convertStringsToChars(precisionRestart);
            self.RestartFile.format         = "\nsampleSize"                                    ...
                                            + "\n" + formatInt                                  ...
                                            + "\nlogSqrtDeterminant"                            ...
                                            + precisionRestart                                  ...
                                            + "\nadaptiveScaleFactorSquared"                    ...
                                            + precisionRestart                                  ...
                                            + "\nmeanVec"                                       ...
                                            + repmat(precisionRestart, 1, self.nd.val)          ...
                                            ...+ "\ncomv_chol"                                  ...
                                            ...+ repmat('\n%.16f', 1, (self.nd.val)^2)          ... do not anymore since comv_covMat is sufficient
                                            ...+ "\ncomv_covMat"                                ...
                                            ...+ repmat(precisionRestart, 1, (self.nd.val)^2)   ...
                                            + "\ncovMat"                                        ...
                                            + repmat(precisionRestart, 1, self.nd.val*(self.nd.val-1)/2 + self.nd.val)      ...
                                            + "\n"                                      ...
                                            ;
            %***********************************************************************************************************************
            
        end % function setupOutputFiles

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end