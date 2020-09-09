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