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
classdef DateTime_class < handle

    properties (Constant)
        CLASS_NAME = "@DateTime_mod"
    end

    properties
        date            = []
        time            = []
        zone            = []
        century         = []
        year            = []
        month           = []
        day             = []
        hour            = []
        minute          = []
        second          = []
        millisecond     = []
        fancyStyleBasic = []
        fancyStyle      = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function query(self)
            dt                      = datetime('now','TimeZone','local','Format','yyyy/MM/dd - HH:mm:ss.SSS Z');
            dt                      = char(dt);
            self.date               = dt([1:4, 6:7, 9:10]);
            self.time               = dt([14:15, 17:18, 20:25]);
            self.zone               = dt(27:31);
            self.century            = dt(1:2);
            self.year               = dt(1:4);
            self.month              = dt(6:7);
            self.day                = dt(9:10);
            self.hour               = dt(14:15);
            self.minute             = dt(17:18);
            self.second             = dt(20:21);
            self.millisecond        = dt(23:25);
            self.fancyStyleBasic    = dt(1:21);
            self.fancyStyle         = strcat(dt, ' UTC');
        end
   
    %***********************************************************************************************************************************
    %*********************************************************************************************************************************** 
        
    end
    
%***********************************************************************************************************************************
%***********************************************************************************************************************************    
    
    methods(Static)
        
    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function getNiceDateTime = getNiceDateTime()
            dt              = datetime('now','TimeZone','local','Format','yyyy/MM/dd HH:mm:ss.SSS Z');
            dt              = char(dt);
            getNiceDateTime = convertCharsToStrings(dt(1:19));
            
            dt              = datetime('now','TimeZone','local','Format','yyyy/MM/dd - HH:mm:ss.SSS Z');
            dt              = char(dt);
            getNiceDateTime = convertCharsToStrings(dt(1:21));
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end