% demoPerturbation.m
%
% Demonstrate usage of the `addEIS` and `simEIS` utilities for running
% perturbation studies in COMSOL.
%
% 2025.02.21 | Update for CDC paper | WH
% 2025.01.01 | Created | GLP

clear; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename))); end
addpath(genpath(fullfile('.','UTILITY')));
addpath(genpath(fullfile('.','XLSX_CELLDEFS')));
addpath(genpath(fullfile('.','TFS')));

SOCs = 0.25:0.25:1;
freqs = logspace(-6,4,200);
CL = lines;
cellData = loadCellParams('cellNMC30.xlsx');

% Generate COMSOL model, add perturbation study for EIS.
model = genFOM(cellData);
model = addEIS(model);

% Perform perturbation study.
[FOMout,model] = simEIS(model,SOCs,freqs);

fprintf('Plotting results...\n');
w = 2*pi*FOMout.freqs; s = 1j*w;
for theZ = 1:length(FOMout.z0)
  figure(50); % Plot COMSOL perturbation impedance as Nyquist
  plot(real(FOMout.runData(theZ).Z),-imag(FOMout.runData(theZ).Z),...
       '-','color',CL(theZ,:)); hold on

  figure(theZ); % Plot COMSOL perturbation impedance as Bode Mag
  semilogx(FOMout.freqs,20*log10(abs(FOMout.runData(theZ).Z))); hold on

  % Get TF impedance
  cellParams = evalSetpoint(cellData,s,FOMout.z0(theZ),298.15);
  if isfield(cellParams.common,'s')
    cellParams.common = rmfield(cellParams.common,'s');
  end    
  [Zcell,Zp] = getImpedance(s,cellParams);
  
  figure(50); % Plot TF impedance as Nyquist
  if theZ == 1, set(gca,'colororderindex',1); end
  plot(real(Zcell),-imag(Zcell),'--'); grid on

  figure(theZ); % Plot TF impedance as Bode
  semilogx(FOMout.freqs,20*log10(abs(Zcell)),'--'); grid on
  title(sprintf('Bode magnitude for SOC = %d%%',round(FOMout.z0(theZ)*100)));
  xlabel('Frequency (Hz)'); ylabel('Mag (dB)')
  thesisFormat;
end

function [Z,Zphase] = getImpedance(s,cellParams)
  [phise_tf,~] = tfPhiseInt(s,[0,3],cellParams);
  [phie_tf,~]  = tfPhie(s,3,cellParams);
  if isfield(cellParams.const,'Rc')
    Rc = cellParams.const.Rc;
  else
    Rc = 0;
  end
  
  Z = -(phise_tf(2,:) - phise_tf(1,:) + phie_tf) + Rc;
  Zphase = angle(Z);
  Zphase(Zphase>0) = Zphase(Zphase>0) - 2*pi;
end