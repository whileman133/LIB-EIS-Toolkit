% function [thetaeTF,aux] = tfThetae(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   thetaeTF = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the TF for normalized concentration of lithium in
% the electrolyte.
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

function [thetaeTF,aux] = tfThetae(s,locs,cellData) %#ok<*NASGU>
  % Initialize outputs
  s = s(:).';          % force "s" to be a row vector
  thetaeTF = zeros(length(locs),length(s));
  aux.hfGain = 0*locs(:);                % high-frequency gain
  aux.dcGain = dcGain(locs);
  aux.res0 = 0*locs(:);                  % integrator residue
  aux.names = cell(length(locs),1);      % TF names
  aux.xLoc = zeros(1,length(locs));      % TF locations  
  
  % Get common components... (next three lines needed to keep MATLAB happy)
  L1n=[];L2n=[];L1s=[];L1p=[];L2p=[];c1n=[];c2n=[];c3n=[];c4n=[];c1s=[]; 
  c2s=[];c1p=[];c2p=[];c3p=[];c4p=[];j1n=[];j2n=[];j3n=[];j4n=[];
  j1p=[];j2p=[];j3p=[];j4p=[];Zsen=[];Zsep=[];Zsn=[];Zsp=[];Isn=[];Isp=[]; 

  [C,L,~,~] = tfCommon(s,cellData); 
  eval(cellData.common.ind); % set c##,L##,J##,Z## indices

  % Compute TF for each cell region
  indNeg = find(locs <= 1);
  if ~isempty(indNeg)
    for k = 1:length(indNeg)
      z = locs(indNeg(k));
      thetaeTF(indNeg(k),:) = ...
          C(c1n,:).*exp(L(L1n,:)*(z-1)) + C(c2n,:).*exp(-L(L1n,:)*z) + ...
          C(c3n,:).*exp(L(L2n,:)*(z-1)) + C(c4n,:).*exp(-L(L2n,:)*z);
      aux.names{indNeg(k)} = 'negThetae';
      aux.xLoc(indNeg(k)) = z;                                                
    end
  end
% stop
  indSep = find(locs>1 & locs<2);
  if ~isempty(indSep)
    for k = 1:length(indSep)
      z = locs(indSep(k))-1;
      thetaeTF(indSep(k),:) = ...
          C(c1s,:).*exp(L(L1s,:)*(z-1)) + C(c2s,:).*exp(-L(L1s,:)*z);
      aux.names{indSep(k)} = 'sepThetae';
      aux.xLoc(indSep(k)) = z+1;                                                        
    end
  end

  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      thetaeTF(indPos(k),:) = ...
          C(c1p,:).*exp(L(L1p,:)*(z-1)) + C(c2p,:).*exp(-L(L1p,:)*z) + ...
          C(c3p,:).*exp(L(L2p,:)*(z-1)) + C(c4p,:).*exp(-L(L2p,:)*z);
      aux.names{indPos(k)} = 'posThetae';
      aux.xLoc(indPos(k)) = 3-z;                                                                
    end
  end
  if ~isempty(find(s==0,1))
    thetaeTF(:,s==0) = aux.dcGain;
  end

  % Compute low-frequency gains (f->0)
  function dcGains = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize

    T = cellData.const.T;

    % Set up constant value
    psi = cellData.const.psi;

    % set up negative-electrode constants...
    kappan = cellData.neg.kappa;
    qen = cellData.neg.qe;

    % set up separator constants
    kappas = cellData.sep.kappa; % not SOC dependent
    qes = cellData.sep.qe;       % not SOC dependent

    % set up positive-electrode constants...
    kappap = cellData.pos.kappa;
    qep = cellData.pos.qe;

    indNeg = find(locs <= 1);
    if ~isempty(indNeg)
      for kk = 1:length(indNeg)
        z = locs(indNeg(kk));
        dcGains(indNeg(kk)) = (qen*kappap*kappas...
            +qep*(2*kappan*kappas+3*kappap*kappas+6*kappan*kappap)...
            +qes*(3*kappap*kappas+3*kappan*kappap))...
            /(6*kappan*kappap*kappas*(psi*T)*(qen+qep+qes))...
            -z^2/(2*kappan*(psi*T));
      end
    end

    indSep = find(locs>1 & locs<2);
    if ~isempty(indSep)
      for kk = 1:length(indSep)
        z = locs(indSep(kk))-1;
        dcGains(indSep(kk)) = (-2*qen*kappap*kappas...
            +qep*(2*kappan*kappas+6*kappan*kappap)+3*qes*kappan*kappap)...
            /(6*kappan*kappap*kappas*(psi*T)*(qen+qep+qes))-z/(kappas*(psi*T));
      end
    end

    indPos = find(locs >= 2);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        dcGains(indPos(kk)) = ...
            (qen*(-2*kappap*kappas-3*kappan*kappas-6*kappan*kappap)...
            -qep*kappan*kappas+qes*(-3*kappan*kappas-3*kappan*kappap))...
            /(6*kappan*kappap*kappas*(psi*T)*(qen+qep+qes))...
            +z^2/(2*kappap*(psi*T));
      end
    end
  end
end