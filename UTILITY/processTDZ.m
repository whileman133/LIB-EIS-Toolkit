function spectra = processTDZ(simData,varargin)
%PROCESSTDZ Process time-domain data collected by EIS simulation in COMSOL.
%
% spectra = PROCESSTDZ(simData) converts the time-domain data collected by
%   COMSOL EIS simulation performed using SIMFOM.m into frequency-domain
%   "transfer-function" quantities. SIMDATA is the output produced by
%   SIMFOM.m. SPECTRA is a structure containing the transfer-function
%   spectra (see description below), which is obtained using
%   the fast-Fourier-transform (FFT) algorithm.
%
% spectra = PROCESSTDZ(...,'NumHarmonics',H) specifies the number of
%   harmonic spectra to compute. H=1 computes the linear spectra only;
%   H=2 computes both linear and second-harmonic spectra; H=3 computes
%   linear, second-, and third-harmonic spectra. Default H=2.
%
% spectra = PROCESSTDZ(...,'EvalLinTF',true) also evalulates the linear
%   transfer-function model for the cell and stores the result in the
%   output for comparison to the linear COMSOL spectra.
%
% spectra = PROCESSTDZ(...,'EvalLinTF',true,'NumTFFreqPoints',N) also
%   specifies the number of frequency-points at which to evalulate the
%   linear transfer-functions. Default N=1000.
%
% The output, SPECTRA, is a structure with the following fields:
%   lin   : structure array of linear spectra for each variable [VariableUnit/A]
%   h2    : (if H>=2) structure array of second-harmonic spectra [VariableUnit/A^2]
%   h3    : (if H>=3) structure array of third-harmonic spectra [VariableUnit/A^3]
%   freq  : frequency vector corresponding to entries in the structure arrays [Hz]
%   xlocs : structure mapping electrochemical variables to the x-locations
%           where the spectra are evalulated (only for internal cell variables)
%   tf    : (present if EvalLinTF==true) structure array of linear spectra
%           evalulated from transfer functions for each variable.
%   tfFreq: (present if EvalLinTF==true) frequency vector corresponding to
%           the tf structure array [Hz].
%   param : structure of parameter values supplied to the function.
%
% -- Examples --
% Let:
%   spectra = processTDZ(simData);
%
% 1. Fetch the linear impedance spectrum:
%    Z1 = [spectra.lin.Zcell];  % row vector, dim1=freq
%
% 2. Fetch the linear Thetae(jw)/Iapp(jw) spectrum at x=3:
%    ind = find(spectra.xlocs.Thetae==3,1,'first')
%    Thetae = [spectra.lin.Thetae]; % matrix, dim1=x-loc, dim2=freq
%    Thetae3 = Thetae(ind,:); % row vector, dim1=freq
%
% -- Changelog --
% 2025.02.21 | Updates for CDC paper | WH
% 2023.04.05 | Created | Wesley Hileman <whileman@uccs.edu>

isinteger = @(x)floor(x)==x;
parser = inputParser;
parser.addRequired('simData',@(x)isstruct(x)&&strcmp(x.origin__,'simTDZ'));
parser.addParameter('NumHarmonics',2,@(x)isscalar(x)&&isinteger(x)&&x>=1);
parser.addParameter('EvalLinTF',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('NumTFFreqPoints',1000,@(x)isscalar(x)&&isinteger(x));
parser.addParameter('Debug',false,@(x)isscalar(x)&&islogical(x));
parser.parse(simData,varargin{:});
p = parser.Results; % structure of validated arguments
H = p.NumHarmonics; % number of harmonics to retain

varNames = setdiff(fieldnames(simData.ss),{'param','time'});
harmNames = arrayfun(@(x)sprintf('h%d',x),1:H,'UniformOutput',false);
harmNames{1} = 'lin';

for h = 1:H
  spectra.(harmNames{h}) = [];
end

for k = length(simData.ss):-1:1
  freqkData = simData.ss(k);

  % Construct FFT frequency vector.
  fs = freqkData.param.fs;  % sampling rate [Sa/s]
  N = freqkData.param.N;    % number of samples available for fft
  f0 = freqkData.param.f0;  % fundemental frequency [Hz]
  tfFreq = 0:(fs/N):(fs/2); % frequency vector [Hz]
  fh = (1:H)*f0;            % harmonic frequencies [Hz]
  [resid,ih] = min(abs(fh-tfFreq')); % harmonic indicies into fft vector
  if any(resid./fh>1e-6)
    warning('Possible FFT leakage for f=%.5g!',f0);
  end

  % Compute frequency spectra of the electrochemical variables.
  Vars = struct;
  if p.Debug
    VarsFullData = struct;
  end
  for j = 1:length(varNames)
    varName = varNames{j};
    timePosData = freqkData.(varName);
    freqPosData = fft(timePosData);
    % Use single-sided spectrum.
    freqPosData = 2*freqPosData(1:length(tfFreq),:)/size(freqPosData,1);
    % Retain only the 1:H harmonics.
    Vars.(varName) = freqPosData(ih,:);
    if p.Debug
      VarsFullData.(varName) = freqPosData;
    end
  end % for

  % Debug output.
  if p.Debug && k==1
    % Time-domain quantities.
    time = freqkData.time;
    iapp = freqkData.Iapp;
    vcell = freqkData.Vcell;
    vcelltilde = vcell-mean(vcell);
    % Frequency-domain quantities.
    Iapp = VarsFullData.Iapp;
    Vcell = VarsFullData.Vcell;

    figure;
    plot(time,iapp,'DisplayName','$i_\mathrm{app}(t)$ [A]'); hold on; 
    plot(time,vcelltilde*100,'DisplayName','$(v_\mathrm{cell}(t) - U_\mathrm{ocv})\times 100$ [V]');
    xlim([min(time) max(time)]);
    ylim([min(iapp) max(iapp)]*1.7);
    legend('Location','best','Interpreter','latex');
    xlabel('Time, $t$ [s]','Interpreter','latex');
    ylabel('Applied Current or Cell Voltage');
    title( ...
      sprintf('Time-Domain Response: $f=%.0f\\,\\mathrm{mHz}$',f0*1000), ...
      'Interpreter','latex');
    thesisFormat;
    figure;
    semilogx(tfFreq,20*log10(abs(Iapp)),'-o','DisplayName','$I_\mathrm{app}(f)$'); hold on;
    semilogx(tfFreq,20*log10(abs(Vcell)),'-o','DisplayName','$V_\mathrm{cell}(f)$');
    xlim([min(tfFreq) max(tfFreq)]);
    xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
    ylabel('Magnitude [dB]');
    title( ...
      sprintf('Frequency Spectra: $f=%.0f\\,\\mathrm{mHz}$',f0*1000), ...
      'Interpreter','latex');
    legend('Location','best','Interpreter','latex');
    thesisFormat('LineMarkerSize',4,'LineLineWidth',1.3);
  end % if

  % Normalize to input current phasor to obtain "transfer function" form.
  II = Vars.Iapp(1).^(1:H).';
  for j = 1:length(varNames)
    varName = varNames{j};
    Vars.(varName) = Vars.(varName)./II;
  end

  % Compute cell impedance.
  Vars.Zcell = -Vars.Vcell;  % Don't forget the negative sign!

  % Store results.
  varsn = fieldnames(Vars);
  for j = 1:length(varsn)
    varName = varsn{j};
    freqPosData = Vars.(varName);
    for h = 1:H
      % Use column vector for variables with more than so than
      % one x-location so that it's easier to collect results using
      % the structure-array notation [structureArray.fieldName].
      spectra.(harmNames{h})(k).(varName) = freqPosData(h,:).';
    end
  end
end

spectra.freq  = simData.param.freq;
if isfield(simData,'xlocs')
  spectra.xlocs = simData.xlocs;
else
  % An earlier version of simEIS did not include the xlocs field.
  spectra.xlocs = simData.param.Vars;
end
spectra.xlocs.Iapp  = [];
spectra.xlocs.Vcell = [];
spectra.xlocs.Zcell = [];

if p.EvalLinTF
  % Compare linear COMSOL spectra to those predicted by transfer functions.

  % Structure mapping COMSOL variables to TF functions.
  % See end of this file for definitions of custom functions!
  var2tf.Thetae = @tfThetae;
  var2tf.Idl = @tfIdl;
  var2tf.If = @tfIf;
  var2tf.Ifdl = @tfIfdl;
  var2tf.Phie = @getPhieTF;    % ground ref. at x=0- (at neg cc)
  var2tf.PhieTilde = @tfPhie;  % ground ref. at x=0+ (in electrolyte)
  var2tf.Phis = @tfPhis;
  var2tf.Phise = @tfPhiseInt;
  var2tf.Thetass = @tfThetassInt;
  var2tf.Vcell = @getVcellTF;
  var2tf.Zcell = @getZcellTF;

  cellModel = simData.param.cellModel;
  socPct = simData.param.socPct;
  TdegC = simData.param.TdegC;
  tfFreq = logspace( ...
    log10(min(spectra.freq)),log10(max(spectra.freq)),p.NumTFFreqPoints);
  ss = 1j*2*pi*tfFreq;
  mod = evalSetpoint(cellModel,ss,socPct/100,TdegC+273.15);

  varNames = setdiff(fieldnames(spectra.lin),{'Iapp'});
  spectra.tf = [];
  for k = 1:length(varNames)
    varName = varNames{k};
    if ~isfield(var2tf,varName)
      continue;
    end
    xlocs = spectra.xlocs.(varName);
    tfValues = var2tf.(varName)(ss,xlocs,mod);
    for j = length(ss):-1:1
      spectra.tf(j).(varName) = tfValues(:,j);
    end
  end

  spectra.tfFreq = tfFreq;
end

spectra.param = p;
spectra.origin__ = 'processTDZ';

end

function PhieTF = getPhieTF(s,x,cellModel)
% Phie with ground reference at x=0- (at negative electrode
% current-collector).
Phise0 = tfPhiseInt(s,0,cellModel);
PhieTilde = tfPhie(s,x,cellModel);
PhieTF = PhieTilde - Phise0;
end

function VcellTF = getVcellTF(s,~,cellModel)
Phise = tfPhiseInt(s,[0 3],cellModel);
PhieTilde3 = tfPhie(s,3,cellModel);
VcellTF = -Phise(1,:) + Phise(2,:) + PhieTilde3 - cellModel.const.Rc;
end

function ZcellTF = getZcellTF(s,~,cellModel)
% Negate cell-voltage TF!
ZcellTF = -getVcellTF(s,[],cellModel);
end
