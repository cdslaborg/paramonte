classdef OS_class < handle

    properties (Constant)
        CLASS_NAME  = "@OS_class"
    end

    properties
        name        = []
        slash       = []
        isWindows   = false
        Err         = Err_class()
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)
    
    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function queryOS(self)
            
            FUNCTION_NAME   = self.CLASS_NAME + "@queryOS()";
            
            self.Err.occurred = false;
            self.Err.msg    = "";

            if ispc
                self.name           = "WINDOWS";
                self.isWindows      = true;
                self.slash = "\";
            elseif ismac
                self.name           = "MAC";
                self.isWindows      = false;
                self.slash          = "/";
            elseif isunix 
                self.name           = "UNIX";
                self.isWindows      = false;
                self.slash          = "/";
            else
                self.name           = "";
                self.isWindows      = false;
                self.Err.occurred   = false;
                self.Err.msg        = FUNCTION_NAME + ": Unknown error occurred while attempting to read "...
                                    + "the Operating System's name";
            end



            
            % [self.name, self.Err]   = EnvVar_class.getEnvVar("OS");
            
            % if self.Err.occurred
                % self.Err.msg    = FUNCTION_NAME + ": Error occurred while querying OS type." + newline + self.Err.msg;
                % self.name       = "";
                % return
            % end
            
            % self.name       = strtrim(self.name);
            
            % if length(self.name) >= 7
                % if lower(self.name(1:7)) == "windows"
                    % self.isWindows  = true;
                    % self.slash      = "\";
                % end
            % else % it's either Linux- or Darwin- based OS
            
                % self.isWindows  = false;
                % self.slash      = "/";
                
                % RFN = RandomFileName_class([],"queryOS",[]);
                % if RFN.Err.occurred
                    % self.Err        = RFN.Err;
                    % self.Err.msg    = FUNCTION_NAME + ": Error occurred while inquiring OS type." + newline + self.Err.msg;
                    % self.name = "";
                    % return
                % end
                
                % [self.Err.occured, self.Err.msg] = system("uname > " + RFN.path);
                % if self.Err.occurred
                    % self.Err.msg    = FUNCTION_NAME + ": Error occurred while executing command 'uname > " + RFN.path + "'." + newline...
                                    % + self.Err.msg;
                    % self.name       = "";
                    % return
                % end
                
                % [fid, self.Err.msg] = fopen("RFN.path","r");
                % if fid > 2
                    % self.Err.occurred = true;
                    % self.Err.msg    = FUNCTION_NAME + ": Error occurred while opening file = '" + RFN.path + "'." + newline...
                                    % + self.Err.msg;
                    % OS.name = "";
                    % return
                % end
                
            % %%%%%%%%%%%     UNFINISHED CODE     %%%%%%%%%%%%
                
            % end
            
            % self.name = strtrim(self.name);
            
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end