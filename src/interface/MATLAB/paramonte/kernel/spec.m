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
classdef spec < handle

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

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Hidden)
    
    % These methods have been implemented to override the default 'handle' class methods, 
    % so that they won't pop-up after pressing 'Tab' button.

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function addlistener    (self)  end
        function delete         (self)  end
        function findobj        (self)  end
        function findprop       (self)  end
        function valid          (self)  end
        function listener       (self)  end
        function notify         (self)  end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end


