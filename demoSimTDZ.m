% demoSimTDZ.m
% 
% Demonstrates the usage of the simTDZ.m utility function for running
% time-domain EIS simulations in COMSOL. Run this file in COMSOL with MATLAB
% LiveLink.
%
% Type `help simTDZ` for detailed usage information.
%
% -- Changelog --
% 2023.04.05 | Created | WH

clear; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename))); end
addpath(genpath(fullfile('.','UTILITY')));
addpath(genpath(fullfile('.','XLSX_CELLDEFS')));
addpath(genpath(fullfile('.','TFS')));

% Constants.
cellFile = 'cellNMC30.xlsx';  % Name of cell parameters spreadsheet.
freq = logspace(-3,5,20);   % Frequency points to evalulate in the spectrum [Hz].
socPct = 50;                % Cell SOC setpoint [%].
TdegC = 25;                 % Cell temperature [degC].
I = 2;                      % Amplitude of Iapp sinusoids [A].

% The following structure specifies which electrochemical variables to
% store in addition to Vcell (field names) as well as the x-locations
% where those variables should be evaluated (field values).
Vars.Phise = [0 3];
Vars.Phie = 3;
Vars.Thetae = 0:3;

% Load cell parameters and run EIS simulation in COMSOL.
cellModel = loadCellParams(cellFile);
simData = simTDZ(cellModel,freq,socPct,TdegC, ...
    'Vars',Vars,'I',I,'Verbose',true);

% Save results to disk.
save("demoTDZ.mat","simData");