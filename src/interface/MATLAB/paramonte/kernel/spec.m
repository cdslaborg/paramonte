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

classdef spec < handle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Hidden)
        url = "https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/";
    end

    properties
        sampleSize                          = []
        randomSeed                          = []
        description                         = []
        outputFileName                      = []
        outputDelimiter                     = []
        chainFileFormat                     = []
        variableNameList                    = []
        restartFileFormat                   = []
        outputColumnWidth                   = []
        overwriteRequested                  = []
        outputRealPrecision                 = []
        silentModeRequested                 = []
        domainLowerLimitVec                 = []
        domainUpperLimitVec                 = []
        parallelizationModel                = []
        progressReportPeriod                = []
        targetAcceptanceRate                = []
        maxNumDomainCheckToWarn             = []
        maxNumDomainCheckToStop             = []
        ...ParaMCMC variable
        chainSize                           = []
        startPointVec                       = []
        sampleRefinementCount               = []
        sampleRefinementMethod              = []
        randomStartPointRequested           = []
        randomStartPointDomainLowerLimitVec = []
        randomStartPointDomainUpperLimitVec = []
        ...ParaDRAM variables
        scaleFactor                         = []
        proposalModel                       = []
        proposalStartCovMat                 = []
        proposalStartCorMat                 = []
        proposalStartStdVec                 = []
        adaptiveUpdateCount                 = []
        adaptiveUpdatePeriod                = []
        greedyAdaptationCount               = []
        delayedRejectionCount               = []
        burninAdaptationMeasure             = []
        delayedRejectionScaleFactorVec      = []
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access=public)

        function helpme(self, specification)
            %   Return help for the input specification. If the input specification is missing,
            %   then help for all specifications will be returned.
            %
            %   Parameters
            %   ----------
            %
            %       specification (optional)
            %           A string or char vector representing the name of the input specification.
            %
            %   Returns
            %   -------
            %
            %       None
            %
            %   Example usage
            %   -------------
            %
            %       pmpd.spec.helpme()              % return help for all specifications.
            %       pmpd.spec.helpme("chainSize")   % return help for the specification chainSize.
            %
            specLower = "";
            errorOccurred = true;
            if nargin==2
                if isstring(specification) || ischar(specification)
                    specList = properties(self);
                    specListLen = length(specList);
                    for i = 1:specListLen
                        if strcmpi(specList{i}, specification)
                            specLower = "#" + lower(specification);
                            errorOccurred = false;
                            break;
                        end
                    end
                end
            elseif nargin==1
                errorOccurred = false;
            end
            if errorOccurred
                error   ( newline ...
                        + "The input specification must be a string or char vector whose value is the name of one of the " ...
                        + "specification properties of the sampler as specified in the ``spec`` component. Example usage: " + newline ...
                        + newline ... 
                        + "    helpme()" + newline ...
                        + "    helpme(""chainSize"")" + newline ...
                        + newline ... 
                        );
            else
                disp("See this page: " + href(self.url+specLower));
                if nargin==1
                    disp( "To get help on a particular simulation specification, try: helpme(specification), like, helpme(""chainSize"")" );
                end
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Hidden)

    % These methods have been implemented to override the default 'handle' class methods, 
    % so that they won't pop-up after pressing 'Tab' button.

        function addlistener    (self)  end
        function delete         (self)  end
        function findobj        (self)  end
        function findprop       (self)  end
        function valid          (self)  end
        function listener       (self)  end
        function notify         (self)  end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


