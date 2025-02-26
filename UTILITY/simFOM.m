% SIMFOM This utility function executes a COMSOL model for specified 
%   input-current and input-temperature profiles, starting at an initial 
%   cell SOC, using the LiveLink for MATLAB interface. The returned COMSOL
%   object can be saved to a file using mphsave.m, or loaded into the 
%   COMSOL GUI using mphlaunch.m. 
%
% Usage 1: [FOM,FOMout] = simFOM(FOM,simData)
% Usage 2: [FOM,FOMout] = simFOM(...,'InputType',typestring)
% Usage 3: [FOM,FOMout] = simFOM(...,'DebugFlag',trueORfalse)
% 
% Inputs:
%   FOM = COMSOL object containing the full-order model, created by 
%     genFOM.m
%   simData = simulation profile loaded using "loadInput" or created
%     programmatically. A structure with the following required fields 
%     depending on the 'InputType' parameter:
%     - 'lut': SOC0, time, T, Iapp
%     - 'sin': SOC0, time, T, freq, mag
%     - 'pwm': SOC0, freq, mag, D, phase, squashTime
%     See below for a detailed description of each field.
%   InputType = type of input to use. typestring should be one of:
%     - 'lut' (default): regular lookup table / interpolation
%     - 'sin': sinusoidal input with specified amplitude and frequency
%     - 'pwm': pulse-width-modulated input
%   DebugFlag = logical value (true or false) indicating whether or not to
%     output status messages to the command window (default true)
%
% Output:
%   FOM     = Updated COMSOL object containing the full-order model and all
%             simulated data  
%   FOMout  = Data structure with results from simulation
%
% simData structure fields:
%   SOC0 = the initial state-of-charge [%]
%   TSHIFT = amount by which to forward-shift input waveform in time [s]
%   time = a vector of simulation time values [s]
%   T = the scalar temperature or vector of temperature values
%     corresponding to the time vector [degC]
%   Iapp = a vector of applied current values corresponding to the
%     time vector [A]
%   freq = the cyclic frequency of the sinusoidal input or the PWM
%     waveform [Hz]
%   mag = the amplitude of the sinusoidal input or PWM waveform [A]
%   D = the duty ratio of the PWM input, 0.5 => 50% duty cycle [-]
%   phase = the phase-offset of the PWM input [rad]
%   squashTime = initial portion of applied current "squashed" [s]
%
% -- Note --
% Electrochemical Impedance Spectroscopy (EIS) Simulations:
%   EIS simulation executes much faster (and more accurately) in COMSOL 
%   when using the 'sin' input type rather than the default 'lut' 
%   (an analytic expression is employed instead of an interpolation 
%   function when using 'sin').
%
% -- Copyright --
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl.
% License, v. 1.0. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L. and Trimboli,
% M. Scott, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2021. 
%
% -- Changelog --
% 2023.04.04 | Add 'InputType' parameter | Wesley Hileman <whileman@uccs.edu>

function [FOM,FOMout] = simFOM(FOM,simData,varargin)
  import com.comsol.model.*      %#ok<NSTIMP>
  import com.comsol.model.util.* %#ok<NSTIMP>
  FOMout = [];

  % Validate input parameters.
  iscomsolmodel = @(x)isa(x,'com.comsol.clientapi.impl.ModelClient');
  issimdata = @(x)isstruct(x);
  isinputtype = @(x)(isstring(x)||ischar(x))&&any(strcmpi(x,{'lut','sin','pwm'}));
  isdebugflag = @(x)isscalar(x)&&islogical(x);
  parser = inputParser;
  parser.addRequired('FOM',iscomsolmodel);
  parser.addRequired('simData',issimdata);
  parser.addParameter('InputType','lut',isinputtype);
  parser.addParameter('DebugFlag',true,isdebugflag);
  parser.parse(FOM,simData,varargin{:});
  p = parser.Results;  % structure of validated parameters
  debugFlag = p.DebugFlag;
  
  % Shift time vector when COMSOL samples all output values by TSHIFT to
  % allow COMSOL to respond to abrupt changes in input current. 
  %
  % If TSHIFT=0, then COMSOL output does not have time to respond to step-
  % like instant changes in input current before sampling voltage (etc.), 
  % so the voltage change does not appear until the following time sample,
  % making the overall voltage profile appear delayed by one time sample.
  % That is, a nonzero TSHIFT is needed to capture feedthrough resistance
  % (ESR) effects. A TSHIFT of 0.1% of the sample period works; a shift of
  % 1% might be slightly better.
  if isfield(simData,'TSHIFT')
    TSHIFT = simData.TSHIFT;
  else
%     TSHIFT = 1e-2; 
    TSHIFT = 0.0025;
  end
                 
  updateInputs;
  runStudy;
  collectResults;
  msg('\n');
  
  function updateInputs
    msg('Updating inputs in FOM...');
    tvec = simData.time(:);
    Tvec = simData.T(:);
    
    % Replace default input-current function
    FOM.func.remove('Ides'); % we'll replace the input function
    if strcmpi(p.InputType,'lut')
        % Lookup table.
        ivec = simData.Iapp(:);
        if length(tvec) == length(ivec),  ivec = ivec(1:end-1); end
        v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = ivec;
        Ides = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
        Ides = eval(sprintf('{ %s }',Ides));
        FOM.func.create('Ides', 'Piecewise');
        FOM.func('Ides').label('Default current profile');
        FOM.func('Ides').set('funcname', 'inputCurrent');
        FOM.func('Ides').set('arg', 't');
        FOM.func('Ides').set('extrap', 'interior');
        FOM.func('Ides').set('smooth', 'contd2');
        FOM.func('Ides').set('smoothzone', '3E-7');
        FOM.func('Ides').set('pieces', Ides);
        FOM.func('Ides').set('argunit', 's');
        FOM.func('Ides').set('fununit', 'A');
    elseif strcmpi(p.InputType,'sin')
        % Sine wave.
        I = simData.mag;
        f0 = simData.freq;
        FOM.func.create('Ides', 'Analytic');
        FOM.func('Ides').label('Default current profile');
        FOM.func('Ides').set('funcname', 'inputCurrent');
        FOM.func('Ides').set('args', 't');
        FOM.func('Ides').set('fununit', 'A');
        FOM.func('Ides').setIndex('argunit', 's', 0);
        FOM.func('Ides').set('expr', sprintf('%g*sin(2*pi*%g*t)',I,f0));
    elseif strcmpi(p.InputType,'pwm')
        tend = tvec(end);
        period = 1/simData.freq;
        D = simData.D;
        phase = simData.phase;
        mag = simData.mag;
        TS = simData.squashTime; % How much delay
        t = TS; VV = [0 TS 0];
        if phase ~= 0
          VV = [VV; TS TS+period/4 mag];
          t = TS + period/4;
        end
        while t <= tend
          VV = [VV; t t+period*D -mag; t+period*D t+period mag];
          t = t+period;
        end
        Ides = sprintf('''%g'' ''%g'' ''%g'';',VV');
        Ides = eval(sprintf('{ %s }',Ides));
        FOM.func.create('Ides', 'Piecewise');
        FOM.func('Ides').label('Default current profile');
        FOM.func('Ides').set('funcname', 'inputCurrent');
        FOM.func('Ides').set('arg', 't');
        FOM.func('Ides').set('extrap', 'interior');
        FOM.func('Ides').set('smooth', 'contd2');
        FOM.func('Ides').set('smoothzone', '3E-7');
        FOM.func('Ides').set('pieces', Ides);
        FOM.func('Ides').set('argunit', 's');
        FOM.func('Ides').set('fununit', 'A');
    end
    
    % Replace default input temperature function
    FOM.func.remove('Temperature');
    if length(Tvec)==1
        % Constant temperature.
        FOM.func.create('Temperature', 'Analytic');
        FOM.func('Temperature').set('funcname', 'inputTemperature');
        FOM.func('Temperature').set('args', 't');
        FOM.func('Temperature').set('argunit', 's');
        FOM.func('Temperature').set('fununit', 'K');
        FOM.func('Temperature').set('expr',sprintf('%g',Tvec+273.15));
    else
        % Time-variant temperature.
        if length(tvec) == length(Tvec),  Tvec = Tvec(1:end-1); end
        v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = Tvec+273.15;
        Temperature = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
        Temperature = eval(sprintf('{ %s }',Temperature));
        FOM.func.create('Temperature', 'Piecewise');
        FOM.func('Temperature').set('funcname', 'inputTemperature');
        FOM.func('Temperature').set('arg', 't');
        FOM.func('Temperature').set('extrap', 'periodic');
        FOM.func('Temperature').set('smooth', 'contd2');
        FOM.func('Temperature').set('smoothzone', '3E-7');
        FOM.func('Temperature').set('pieces', Temperature);
        FOM.func('Temperature').set('argunit', 's');
        FOM.func('Temperature').set('fununit', 'K');
    end

    % Replace default initial SOC
    FOM.param.set('z0', num2str(simData.SOC0/100), 'Initial cell SOC');
    
    % Replace default study values
    tvar = tvec+TSHIFT; 
    FOM.study('std1').feature('time').set('tlist', tvar);      
    FOM.sol('sol1').feature('t1').set('tlist', tvar);
  end   % updateInputs
  function runStudy 
    import com.comsol.model.* 
    import com.comsol.model.util.*    
    msg('\nRunning study...');
    if ~ismac, ModelUtil.showProgress(true); end    
    % mphlaunch(mphtags(FOM));    
    try
      FOM.sol('sol1').runAll;
    catch
      warning('COMSOL did not converge! Results computed will be saved.')
    end    
  end      % runStudy
  function collectResults
    msg('\nCollecting results...');    
    % In negative electrode and at its current collector
    vars = {'phi_s','phi_e','theta_e','if','idl','thetass','ifdl'};
    data_neg = mpheval(FOM,vars,'Edim',1,'Selection',1);
    locs_neg = data_neg.p(:);
    FOMout.negIfdl    = data_neg.d7;
    FOMout.negIf      = data_neg.d4;
    FOMout.negIdl     = data_neg.d5;
    FOMout.negPhis    = data_neg.d1;
    FOMout.negPhise   = data_neg.d1 - data_neg.d2;
    FOMout.negThetass = data_neg.d6;

    FOMout.negIfdl0    = FOMout.negIfdl(:,1); 
    FOMout.negIf0      = FOMout.negIf(:,1);
    FOMout.negIdl0     = FOMout.negIdl(:,1);
    FOMout.negPhise0   = FOMout.negPhise(:,1);
    FOMout.negThetass0 = FOMout.negThetass(:,1);
        
    % In the electrolyte
    data_elect = mpheval(FOM,{'theta_e','phi_e'},'Edim',1);
    locs_elect = data_elect.p(:);
    FOMout.Phie   = data_elect.d2;
    FOMout.Thetae = data_elect.d1;
    FOMout.xLocs.Phie   = locs_elect;
    FOMout.xLocs.Thetae = locs_elect;
    
    % In positive electrode and at its current collector
    data_pos = mpheval(FOM,vars,'Edim',1,'Selection',3);
    locs_pos = data_pos.p(:);
    FOMout.posIfdl    = data_pos.d7;
    FOMout.posIf      = data_pos.d4;
    FOMout.posIdl     = data_pos.d5;
    FOMout.posPhis    = data_pos.d1;
    FOMout.posPhise   = data_pos.d1 - data_pos.d2;
    FOMout.posThetass = data_pos.d6;

    FOMout.posIfdl3    = FOMout.posIfdl(:,end); 
    FOMout.posIf3      = FOMout.posIf(:,end); % GLP: 01/04/23
    FOMout.posIdl3     = FOMout.posIdl(:,end);
    FOMout.posPhise3   = FOMout.posPhise(:,end);
    FOMout.posThetass3 = FOMout.posThetass(:,end);

    FOMout.xLocs.Ifdl    = [locs_neg; locs_pos];
    FOMout.xLocs.If      = [locs_neg; locs_pos];
    FOMout.xLocs.Idl     = [locs_neg; locs_pos];
    FOMout.xLocs.Phis    = [locs_neg; locs_pos];
    FOMout.xLocs.Phise   = [locs_neg; locs_pos];
    FOMout.xLocs.Thetass = [locs_neg; locs_pos];
                          
    % Simulation control and primary output
    temp = mpheval(FOM,'T');
    FOMout.T = temp.d1(:,end)-273.15; % Degrees Celsius
    tvar = mpheval(FOM,'t');  tvar = tvar.d1(:,1);
    FOMout.time = tvar-TSHIFT;
    input = mpheval(FOM,'Ides');   
    FOMout.Ides = input.d1(:,end);
    input = mpheval(FOM,'Iapp');   
    FOMout.Iapp = input.d1(:,end);
    input = mpheval(FOM,'Vcell','Edim',0);
    FOMout.Vcell = input.d1(:,end);
    
    % Other cell and electrode quantities
    % ROMout.cellSOC = zeros(duration,1); 
    data_neg = mpheval(FOM,'thetasavg_neg');
    FOMout.negSOC = data_neg.d1;
    data_pos = mpheval(FOM,'thetasavg_pos');
    FOMout.posSOC = data_pos.d1;
    
    % Save each variable
    FOMout.Ifdl    = [FOMout.negIfdl FOMout.posIfdl];
    FOMout.If      = [FOMout.negIf FOMout.posIf];
    FOMout.Idl     = [FOMout.negIdl FOMout.posIdl];
    FOMout.Phis    = [FOMout.negPhis FOMout.posPhis];
    FOMout.Phise   = [FOMout.negPhise FOMout.posPhise];
    FOMout.Thetass = [FOMout.negThetass FOMout.posThetass];
  end % collectResults
  function msg(theText)
    if debugFlag
      fprintf(theText);
    end
  end % msg
end