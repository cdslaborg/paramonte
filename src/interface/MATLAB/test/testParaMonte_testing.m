%**********************************************************************************************************************************
%**********************************************************************************************************************************
%
%  ParaMonte: plain powerful parallel Monte Carlo library.
%
%  Copyright (C) 2012-present, The Computational Data Science Lab
%
%  This file is part of ParaMonte library. 
%
%  ParaMonte is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, version 3 of the License.
%
%  ParaMonte is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
%
%**********************************************************************************************************************************
%**********************************************************************************************************************************

% clear MATLAB space

clc;
clear all;
close all;
clear classes;
format compact; format long;

% set path to the ParaMonte library

pmlibRootDir = './'; % set this path to the paramonte library root dir
addpath(genpath(pmlibRootDir));

% change MATLAB's working directory to the folder containing this script

filePath = mfilename('fullpath');
[currentDir,fileName,fileExt] = fileparts(filePath); cd(currentDir); % Change working directory to source code directory.
cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.

% define the objective function

NDIM = 4; % number of dimensions of the distribution

function logFunc = getLogFunc(point)
    mean    =   [ 0.0,0.0,0.0,0.0];     % mean of the Multivariate Normal distribution
    covmat  =   [ 1.0,0.5,0.5,0.5 ;     ... covariance matrix of the Multivariate Normal distribution
                , 0.5,1.0,0.5,0.5 ;     ...
                , 0.5,0.5,1.0,0.5 ;     ...
                , 0.5,0.5,0.5,1.0 ]
    logFunc = log(mvnpdf(point,mean,covmat));
end

% create a ParaMonte object

pm = paramonte(); % use paramonte("matlab") to use the MATLAB kernel instead of interface

% create a ParaDRAM simulation object

pmpd = pm.ParaDRAM();

% specifity path to the input specification file. 
% This is optional: all simulation specifications can be 
% also set as attributes of the pmpd.spec component of the object.
% KEEP IN MIND: if you set pmpd.inputFile to any non-empty value, then
% the inputFile will override any values specified via pmpd.spec properties.

pmpd.inputFile = string(currentDir) + "/paramonte.in";

pmpd.runSampler ( NDIM          ... number of dimensions of the objective function
                , getLogFunc    ... the objective function: multivariate normal distribution
                );
