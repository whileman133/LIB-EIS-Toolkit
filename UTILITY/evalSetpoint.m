% function cellParams = evalSetpoint(cellParams,s,soc,T)
% 
% Inputs:
%   cellParams = data structure containing cell model parameter functions,
%                most likely loaded from an Excel spreadsheet
%   s          = vector of 1j*w where w is a vector of frequencies at which 
%                transfer functions will be evaluated. This argument can be
%                empty if these frequencies are unknown, but transfer-
%                function computations are sped up if a frequency vector is
%                specified at this time
%   soc        = cell state-of-charge (between 0 and 1) for the setpoint at
%                which the model is evaluated
%   T          = cell temperature (in K) for the setpoint at which the
%                model is evaluated
%
% Output:
%   cellParams = updated data structure containing original data as well as
%                specific values for parameters at the specified setpoint
%
% This utility function evaluates the cell-model parameter values for a
% specific state-of-charge and temperature setpoint. These parameter values
% are created by evaluating the functions in cellParams.function and are
% stored in cellParams.const, cellParams.neg, cellParams,sep, and
% cellParams.pos. If "s" is not empty, cellParams.common is created to
% store some computations common to evaluating all transfer functions -- if
% this operation is done once at this point it speeds up later transfer-
% function evaluations.
% 
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl.
% License, v. 1.0. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L. and Trimboli,
% M. Scott, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2021. 

% Takes the cellParams.function.[ fields ] functions and evaluates at the
% specific given SOC and T setpoint to produce constant values 
% ("s" is a vector of frequencies to input to tfCommon)
function cellParams = evalSetpoint(cellParams,s,soc,T)
  F = fields(cellParams.function.const);
  for k = 1:length(F)
    fn = cellParams.function.const.(F{k});
    val = fn(soc,T);
    cellParams.const.(F{k}) = val;
  end
  
  F = fields(cellParams.function.neg);
  negSOC = cellParams.function.neg.soc(soc,T);
  for k = 1:length(F)
    fn = cellParams.function.neg.(F{k});
    val = fn(negSOC,T);
    cellParams.neg.(F{k}) = val;
  end
  cellParams.neg.soc = negSOC;
  
  F = fields(cellParams.function.sep);
  for k = 1:length(F)
    fn = cellParams.function.sep.(F{k});
    val = fn(soc,T);
    cellParams.sep.(F{k}) = val;
  end
  
  F = fields(cellParams.function.pos);
  posSOC = cellParams.function.pos.soc(soc,T);
  for k = 1:length(F)
    fn = cellParams.function.pos.(F{k});
    val = fn(posSOC,T);
    cellParams.pos.(F{k}) = val;
  end
  cellParams.pos.soc = posSOC;
  
  cellParams.const.soc = soc;
  cellParams.const.T = T;
  
  % SOC-dependent Ds
  if isfield(cellParams.neg,'Dsref')
    cellParams.neg.Ds = cellParams.neg.Dsref*cellParams.const.F/cellParams.const.R/T...
        *cellParams.neg.soc*(cellParams.neg.soc-1)*cellParams.neg.dUocp;
  end
  if isfield(cellParams.pos,'Dsref')
    cellParams.pos.Ds = cellParams.pos.Dsref*cellParams.const.F/cellParams.const.R/T...
        *cellParams.pos.soc*(cellParams.pos.soc-1)*cellParams.pos.dUocp;
  end
  
  if ~isempty(s)
    cellParams.common = []; % 20200629 force "tfCommon" to overwrite 
    [C,L,J,Z,Rct] = tfCommon(s,cellParams);
    cellParams.common.s = s;
    cellParams.common.C = C;
    cellParams.common.L = L;
    cellParams.common.J = J;
    cellParams.common.Z = Z;
    cellParams.common.Rct = Rct;
    cellParams.common.ind = ['L1n=1;L2n=2;L1s=3;L1p=4;L2p=5;' ...
      'c1n=1;c2n=2;c3n=3;c4n=4;c1s=5;c2s=6;c1p=7;c2p=8;c3p=9;c4p=10;'...
      'j1n=1;j2n=2;j3n=3;j4n=4;j1p=5;j2p=6;j3p=7;j4p=8;' ...
      'Zsen=1;Zsep=2;Zsn=3;Zsp=4;Isn=5;Isp=6;'];
  end
end