% demoProcessTDZ.m
%
% Demonstrates usage of the processEIS.m utility function for converting
% time-domain COMSOL simulation data into frequency spectra.
%
% Type `help processTDZ` for detailed usage information.
%
% -- Changelog --
% 2023.04.05 | Created | WH

clear; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename))); end
addpath(genpath(fullfile('.','UTILITY')));
addpath(genpath(fullfile('.','XLSX_CELLDEFS')));
addpath(genpath(fullfile('.','TFS')));
load('demoTDZ.mat');

% Compute linear spectra from COMSOL data.
% Also evalulate linear TFs of same variables at same x-locs for 
% comparison to COMSOL simulation.
spectra = processTDZ(simData, ...
    'NumHarmonics',1, ...
    'EvalLinTF',true, ...
    'NumTFFreqPoints',1000, ...
    'Debug',true);
Z1sim = [spectra.lin.Zcell];
Z1tf = [spectra.tf.Zcell];
ThetaeSim = [spectra.lin.Thetae];
ThetaeTF = [spectra.tf.Thetae];
xxThetae = spectra.xlocs.Thetae;

% Plot linear impedance (Nyqiust).
figure;
plot(real(Z1tf),-imag(Z1tf)); hold on;
plot(real(Z1sim),-imag(Z1sim),'d');
legend('TF','COMSOL','Location','northwest');
thesisFormat;

% Plot linear impedance (Bode).
figure;
loglog(spectra.tfFreq,abs(Z1tf)); hold on;
loglog(spectra.freq,abs(Z1sim),'d');
xlabel('Cyclic Frequency [Hz]');
ylabel('|Z_{cell}| [\Omega]');
legend('TF','COMSOL');
thesisFormat;
figure;
semilogx(spectra.tfFreq,angle(Z1tf)*180/pi); hold on;
semilogx(spectra.freq,angle(Z1sim)*180/pi,'d');
xlabel('Cyclic Frequency [Hz]');
ylabel('\angleZ_{cell} [deg]');
legend('TF','COMSOL','Location','northwest');
thesisFormat;

% Plot linear Thetae/Iapp TF (Nyqiust).
labels1 = arrayfun(@(x)sprintf('x=%.0f TF',x),xxThetae, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.0f FOM',x),xxThetae, ...
    'UniformOutput',false);
colors = cool(length(xxThetae));
figure;
for k = 1:length(xxThetae)
    plot(real(ThetaeTF(k,:)),-imag(ThetaeTF(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxThetae)
    plot(real(ThetaeSim(k,:)),-imag(ThetaeSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Theta_e(j\omega)/Iapp(j\omega)');
xlabel('Real');
ylabel('-Imag');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
thesisFormat([0.1 0.1 0.1 0.1]);