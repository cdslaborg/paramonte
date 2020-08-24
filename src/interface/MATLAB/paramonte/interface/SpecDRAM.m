%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it 
%   under the terms of the GNU Lesser General Public License as published 
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see, 
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms, 
%   if you use any parts of this library for any purposes, 
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
classdef SpecDRAM < SpecMCMC

    properties(Hidden)
        url
    end

    properties
        adaptiveUpdateCount                  = [];
        adaptiveUpdatePeriod                 = [];
        greedyAdaptationCount                = [];
        delayedRejectionCount                = [];
        burninAdaptationMeasure              = [];
        delayedRejectionScaleFactorVec       = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = SpecDRAM(domain)
            self.url = domain + "notes/usage/paradram/specifications/";
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                        + "The input specification must be a string or char vector whose value is the name of one of the properties of the SpecDRAM class. Example usage: " ...
                        + newline + newline ... 
                        + "    helpme()" + newline ...
                        + "    helpme(""chainSize"")" + newline ...
                        + newline + newline ... 
                        );
            else
                disp("See this page: " + href(self.url+specLower));
                if nargin==1
                    disp( "To get help on a particular simulation specification, try: helpme(specification), like, helpme(""chainSize"")" );
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % classdef SpecDRAM
