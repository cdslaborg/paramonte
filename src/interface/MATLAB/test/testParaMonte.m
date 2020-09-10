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

%C:\Users\<username>\AppData\Roaming\MathWorks\MATLAB
%cd D:\Dropbox\Projects\20180101_ParaMonte\git\src\interface\MATLAB\bin
%mex -nojvm CC=icl paramonte.c libparamonte_dynamic_heap_testing_intel_c_windows_x64_mt.lib -output libparamonte_dynamic_heap_testing_intel_windows_x64_mt
%[status,cmdout] = system('matlab -nosplash -nojvm -r "testParaMonte,exit"','-echo');
%clc;
clear all;
close all;
fclose('all');
clear classes;
format compact; format long;
pmlibRootDir = './'; % set this path to the paramonte library root dir
addpath( genpath(pmlibRootDir) ) 
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath); % Change working directory to source code directory.
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.

%-----------------------------------------------------------------------------------------------------------------------------------------------------------
%pm = paramonte("matlab");
pm = paramonte(); % "matlab"
pmpd = pm.ParaDRAM();
%pmpd = ParaDRAM_class();
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
% return
% %pmpd.inputFile = './paramonte.in';

% pmpd.spec.chainSize = 20000;
% %pmpd.spec.adaptiveUpdateCount = 2100000000;
% %pmpd.spec.adaptiveUpdatePeriod = 3000;
% %pmpd.spec.startPointVec = -10;
% pmpd.spec.randomSeed = 35671;
% pmpd.spec.proposalModel = "normal";
% pmpd.spec.targetAcceptanceRate = 0.2;
% pmpd.runSampler(2,@getLogFunc); %@(point)-0.5*sum(point.^2)); %
% %pmpd.spec.outputFileName = "D:\Dropbox\Projects\20180101_ParaMonte\git\src\interface\MATLAB\test\ParaDRAM_run_300420_012424_780";
% %pmpd.readMarkovChain();
% pmpd.readChain(); %"ParaDRAM_run_180420_151344_607";%"D:\Dropbox\Projects\20180101_ParaMonte\git\src\interface\MATLAB\test\ParaDRAM_run_160420_023054_530");
% %pmpd.readSample();
% return
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaMonte variables...
%pmpd.spec.sampleSize                            = 100;                         % Works
%pmpd.spec.randomSeed                            = 7;                            % Works
pmpd.spec.description                           = "Hi there";                   % Works
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
flag = 0;
flag = 1;
if flag
    file    = "./out/temp/A_A_A_A";
    if exist(file + "_process_1_chain.txt"),     delete(file + "_process_1_chain.txt");     end
    if exist(file + "_process_1_progress.txt"),  delete(file + "_process_1_progress.txt");  end
    if exist(file + "_process_1_report.txt"),    delete(file + "_process_1_report.txt");    end
    if exist(file + "_process_1_restart.txt"),   delete(file + "_process_1_restart.txt");   end
    if exist(file + "_process_1_sample.txt"),    delete(file + "_process_1_sample.txt");    end
    pmpd.spec.outputFileName                    = file;
else
    pmpd.spec.outputFileName                    = "./out/temp/";
end
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
%pmpd.spec.outputDelimiter                       = "|";                          % Works
%pmpd.spec.chainFileFormat                       = "compact";                    % Works
%pmpd.spec.variableNameList                      = ["Variable-X", "Variable-Y"]; % Works
%pmpd.spec.restartFileFormat                     = "ASCII";                      % Works
%pmpd.spec.outputColumnWidth                     = 25;                           % Works
%pmpd.spec.outputRealPrecision                   = 17;                           % Works
%pmpd.spec.silentModeRequested                   = 0;                            % Works
%pmpd.spec.domainLowerLimitVec                   = -1.e300 * ones(5,1);%[-4,-4];                      % Works
%pmpd.spec.domainUpperLimitVec                   = +1.e300 * ones(5,1);%[+4,+4];                      % Works
%pmpd.spec.progressReportPeriod                  = 100;                          % Works
%pmpd.spec.targetAcceptanceRate                  = 0.5;                          % Works
% pmpd.spec.maxNumDomainCheckToWarn               = [];
% pmpd.spec.maxNumDomainCheckToStop               = [];
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaMCMC variables...
pmpd.spec.chainSize                             = 1000;                         % Works
%pmpd.spec.startPointVec                         = [0.6,1.2];                    % Works
%pmpd.spec.sampleRefinementCount                 = 1;                            % Works
%pmpd.spec.sampleRefinementMethod                = "someRandomName";             % Works
%pmpd.spec.randomStartPointRequested             = 1;                            % Works
%pmpd.spec.randomStartPointDomainLowerLimitVec   = [0.5, 1.0];                   % Works
%pmpd.spec.randomStartPointDomainUpperLimitVec   = [1.0, 1.5];                   % Works
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    ...ParaDRAM variables...
%pmpd.spec.scaleFactor                           = "1.2 * gelman";               % Works
%pmpd.spec.proposalModel                         = "normal";                     % Works
%pmpd.spec.proposalStartCovMat                   = [0.5, 0.2; 0.1, 0.3];         % Works
%pmpd.spec.proposalStartCorMat                   = [0.4, 0.1; 0.2, 0.3];         % Works
%pmpd.spec.proposalStartStdVec                   = [1, 1];                       % Works
%pmpd.spec.adaptiveUpdateCount                   = 2;                            % Works
%pmpd.spec.adaptiveUpdatePeriod                  = 25;                           % Works
%pmpd.spec.greedyAdaptationCount                 = 2;                            % Works
%pmpd.spec.burninAdaptationMeasure               = 0.5;                          % Works
%pmpd.spec.delayedRejectionCount                 = 2;                            % Works
%pmpd.spec.delayedRejectionScaleFactorVec        = [3, 4];                       % Works
%pmpd.spec.delayedRejectionCount             = 10;
%pmpd.spec.delayedRejectionScaleFactorVec    = [1.1,2,3,1.5,2.1,3,2,1.1,1.2,1.3];
%-----------------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ndim = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return
pmpd.runSampler(ndim, @getLogFunc);
fclose('all');

% scatter(pmpd.Chain.State(1,:),pmpd.Chain.State(2,:));
% plot3(pmpd.Chain.State(1,:),pmpd.Chain.State(2,:),pmpd.Chain.LogFunc);

% for i = 1 : ndim
% 
%     figure;
%     histfit(pmpd.Chain.State(i,:));
%     
%     figure;
%     plot(pmpd.Chain.State(i,:));
%     set(gca,'xscale','log')
%   
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%system(pmpd.LogFile.Path.modified);
%system(pmpd.TimeFile.Path.modified);
%system(pmpd.ChainFile.Path.modified);
%system(pmpd.SampleFile.Path.modified);
%system(pmpd.RestartFile.Path.modified);