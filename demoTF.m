% demoTF.m
%
% Demonstrates usage of the evalSetpoint() and tfXX() functions for
% computing the frequency response of an LIB cell.
%
% -- Changelog --
% 2025.02.25 | Created | WH

clear; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename))); end
addpath(genpath(fullfile('.','UTILITY')));
addpath(genpath(fullfile('.','XLSX_CELLDEFS')));
addpath(genpath(fullfile('.','TFS')));

% Define constants.
cellParams = loadCellParams('cellNMC30.xlsx');
socPct = 50;  % cell SOC setpoint
TdegC = 25;   % cell temperature setpoint
freq = logspace(-3,5,50);  % frequency vector [Hz]
tfs = [
  struct('name','Thetae','latex','{\tilde{\Theta}_\mathrm{e}}','unit','\mathrm{A}^{-1}','fn',@tfThetae,'xlocs',0:3)
  struct('name','Idl','latex','{I_\mathrm{dl}}','unit','\mathrm{A}\,\mathrm{A}^{-1}','fn',@tfIdl,'xlocs',0:3)
  struct('name','Ifdl','latex','{I_\mathrm{fdl}}','unit','\mathrm{A}\,\mathrm{A}^{-1}','fn',@tfIfdl,'xlocs',0:3)
  struct('name','Phie','latex','{\tilde{\Phi}_\mathrm{e}}','unit','\mathrm{V}\,\mathrm{A}^{-1}','fn',@tfPhie,'xlocs',1:3)
  struct('name','Phis','latex','{\tilde{\Phi}_\mathrm{s}}','unit','\mathrm{V}\,\mathrm{A}^{-1}','fn',@tfPhis,'xlocs',0:3)
  struct('name','Phise','latex','{\tilde{\Phi}_\mathrm{se}}','unit','\mathrm{V}\,\mathrm{A}^{-1}','fn',@tfPhiseInt,'xlocs',0:3)
  struct('name','Thetass','latex','{\tilde{\Theta}_\mathrm{ss}}','unit','\mathrm{A}^{-1}','fn',@tfThetassInt,'xlocs',0:3)
  struct('name','Zcell','latex','Z_\mathrm{cell}','unit','\Omega','fn',@getZcellTF,'xlocs',[])
];  % struct array of TFs to evalulate
tfs(strcmp({tfs.name},'Phise')).inset = {[1.9e-4,2.1e-4],[-2e-4,-8e-4],'YSpan',[-1e-5 3e-5]};
plotColors = colorblind10;
plotOpts = {'LineLineWidth',1.5,'LineMarkerSize',4};
plotSave = false;

% Get parameter values at specified SOC / temperature setpoint.
ss = 1j*2*pi*freq;
paramsSP = evalSetpoint(cellParams,ss,socPct/100,TdegC+273.15);

% Evalulate TFs over frequency vector.
for k = 1:length(tfs)
    tf = tfs(k);
    xlocs = tf.xlocs;
    tfValues = tf.fn(ss,xlocs,paramsSP);
    tfs(k).freq = freq;
    tfs(k).values = tfValues.';
end % for

% -- Plotting --
for k = 1:length(tfs)
  tf = tfs(k);
  labels = arrayfun( ...
    @(x)sprintf('$\\tilde{x}=%.2f$',x),tf.xlocs,'UniformOutput',false);
  figure; 
  colororder(plotColors);
  plot(real(tf.values),-imag(tf.values),'o-');
  if isempty(tf.xlocs)
    title( ...
      sprintf('Nyquist: $%s(j\\omega)/I_\\mathrm{app}(j\\omega)$',tf.latex), ...
      'Interpreter','latex');
    xlabel( ...
      sprintf("$%s'(j\\omega)$ [$%s$]",tf.latex,tf.unit), ...
      'Interpreter','latex');
    ylabel( ...
      sprintf("$-%s''(j\\omega)$ [$%s$]",tf.latex,tf.unit), ...
      'Interpreter','latex');
  else
    title( ...
      sprintf('Nyquist: $%s(\\tilde{x},j\\omega)/I_\\mathrm{app}(j\\omega)$',tf.latex), ...
      'Interpreter','latex');
    xlabel( ...
      sprintf("$G_%s'(\\tilde{x},j\\omega)$ [$%s$]",tf.latex,tf.unit), ...
      'Interpreter','latex');
    ylabel( ...
      sprintf("$-G_%s''(\\tilde{x},j\\omega)$ [$%s$]",tf.latex,tf.unit), ...
      'Interpreter','latex');
  end
  setAxesNyquist;
  thesisFormat(plotOpts{:});
  if isfield(tf,'inset') && ~isempty(tf.inset)
    ax = addInset(tf.inset{:});
    setAxesNyquist('axes',ax,'padxPct',0,'padyPct',0);
  end
  if ~isempty(labels)
    legend(labels{:},'Location','best','Interpreter','latex');
  end
  if plotSave
    print(fullfile('plots','tfs',tf.name),'-depsc');
    print(fullfile('plots','tfs',tf.name),'-dpng');
  end
end % for

% Utility function for computing cell impedance.
function ZcellTF = getZcellTF(s,~,cellModel)
    Phise = tfPhiseInt(s,[0 3],cellModel);
    PhieTilde3 = tfPhie(s,3,cellModel);
    ZcellTF = cellModel.const.Rc - (Phise(2,:) + PhieTilde3 - Phise(1,:));
end