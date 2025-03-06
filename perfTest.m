% perfTest.m
% 
% Measure the wall-clock time required to perform impedance calculations.
%
% -- Changelog --
% 2025.03.06 | Created | WH

clear; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename))); end
addpath(genpath(fullfile('.','UTILITY')));
addpath(genpath(fullfile('.','XLSX_CELLDEFS')));
addpath(genpath(fullfile('.','TFS')));

% Constants.
cellParams = loadCellParams('cellNMC30.xlsx');
freq = logspace(-3,4,10);   % Frequency points to evalulate in the spectrum [Hz].
socPct = 50;                % Cell SOC setpoint [%].
TdegC = 25;                 % Cell temperature [degC].
I = cellParams.function.const.Q()/10; % Amplitude of Iapp sinusoids [A].
numTrials = 10;

% Time domain simulation.
wctTDZ = zeros(1,numTrials);
for t = 1:numTrials
  tic;
  simData = simTDZ(cellParams,freq,socPct,TdegC, ...
      'Vars',struct, ...  % don't collect any variables except Vcell
      'I',I, ...
      'Verbose',true);
  spectra = processTDZ(simData, ...
      'NumHarmonics',1, ...
      'EvalLinTF',false, ...
      'Debug',false);
  wctTDZ(t) = toc;
end

% Perturbation.
wctPert = zeros(1,numTrials);
for t = 1:numTrials
  tic;
  model = genFOM(cellParams);
  model = addEIS(model);
  [FOMout,model] = simEIS(model,SOCs,freqs);
  wctPert(t) = toc;
end

% TF.
wctTF = zeros(1,numTrials);
ss = 1j*2*pi*freq;
for t = 1:numTrials
  tic;
  paramsSP = evalSetpoint(cellParams,ss,socPct/100,TdegC+273.15);
  Zcell = getZcellTF(ss,[],paramsSP);
  wctTF(t) = toc;
end

% Save results to disk.
save("perfTest.mat", ...
  "wctTDZ","wctPert","wctTF", ...
  "cellParams","freq","socPct","TdegC","I","numTrials");


% Utility function for computing cell impedance.
function ZcellTF = getZcellTF(s,~,cellModel)
    Phise = tfPhiseInt(s,[0 3],cellModel);
    PhieTilde3 = tfPhie(s,3,cellModel);
    ZcellTF = cellModel.const.Rc - (Phise(2,:) + PhieTilde3 - Phise(1,:));
end