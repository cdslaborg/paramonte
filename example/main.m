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
%   we ask you to acknowledge the ParaMonte library's usage
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear MATLAB space

clc;
clear all;
close all;
clear classes;
format compact; format long;

% set path to the ParaMonte library

pmlibRootDir = './'; % change this path to the ParaMonte library root directory
addpath(genpath(pmlibRootDir));

% change MATLAB's working directory to the folder containing this script

filePath = mfilename('fullpath');
[currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir);
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.

% instantiate an object of class logfunc containing:
%     - the objective function (getLogFunc)
%     - its number of dimensions (NDIM)

logFunc = logfunc(); 

% create a ParaMonte object:
% NOTE: MatDRAM is a pure-MATLAB implementation of the ParaMonte library's ParaDRAM routine. 
% NOTE: To use the MatDRAM kernel instead of the ParaMonte MATLAB interface, try:
% pm = paramonte("matdram");
% instead of the following:
pm = paramonte();

% create a ParaDRAM simulation object

pmpd = pm.ParaDRAM();

% specify the path to the input specification file.
% This is optional: all simulation specifications can be 
% also set as attributes of the pmpd.spec component of the object.
% KEEP IN MIND: if you set pmpd.inputFile to any non-empty value, then
% the inputFile will override any values specified via pmpd.spec properties.
% uncomment the following line to specify the simulation input specifications 
% solely from within the external input file:
% pmpd.inputFile = string(currentDir) + "/paramonte.in";

% NOTE: The following specifications will be ignored if the input file above is uncommented.
% NOTE: See the following link for the complete list of simulation specifications:
% NOTE:
% NOTE:     https://www.cdslab.org/paramonte/notes/usage/paradram/specifications/

pmpd.spec.chainSize = 30000;            % the number of uniquely-sampled points 
pmpd.spec.outputFileName = "./out/";    % only the output folder specified here in the above, implying that 
                                        % the filenames are to be generated automatically by the sampler.

% run the ParaDRAM simulation

pmpd.runSampler ( logFunc.NDIM  ... number of dimensions of the objective function
                , @logFunc.get  ... the objective function: multivariate normal distribution
                );
