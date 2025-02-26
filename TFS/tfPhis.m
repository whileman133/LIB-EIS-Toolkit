% function [phisTF,aux] = tfPhis(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   phisTF   = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the solid-potential TF.
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

function [phisTF,aux] = tfPhis(s,locs,cellData)
  % Set up some constants
  F = cellData.const.F; %#ok<*NASGU>
  T = cellData.const.T;

  % Initialize shared variables needed for calculating gains
  sigman = cellData.neg.sigma;
  sigmap = cellData.pos.sigma;
  
  % Initialize outputs
  s = s(:).';          % force "s" to be a row vector
  [~,L,J,~,Rct] = tfCommon(s,cellData); 

  phisTF = zeros(length(locs),length(s));
  aux.hfGain = calcDterm(locs);          % high-frequency gain
  aux.dcGain = dcGain(locs);             % low-frequency gain
  aux.res0 = 0*locs(:);                  % integrator residue
  aux.names = cell(length(locs),1);      % TF names
  aux.xLoc = zeros(1,length(locs));      % TF locations

  % Get common components... (next three lines needed to keep MATLAB happy)
  L1n=[];L2n=[];L1s=[];L1p=[];L2p=[];c1n=[];c2n=[];c3n=[];c4n=[];c1s=[]; 
  c2s=[];c1p=[];c2p=[];c3p=[];c4p=[];j1n=[];j2n=[];j3n=[];j4n=[];
  j1p=[];j2p=[];j3p=[];j4p=[];Zsen=[];Zsep=[];Zsn=[];Zsp=[];Isn=[];Isp=[];
  eval(cellData.common.ind); % set c##,L##,J##,Z## indices

  % Compute TF at each electrode
  indNeg = find(locs <= 1);
  if ~isempty(indNeg)
    for k = 1:length(indNeg)
      z = locs(indNeg(k));
      phisTF(indNeg(k),:) = (-z ...
          + J(j1n,:).*(exp( L(L1n,:)*(z-1))-exp(-L(L1n,:)) ...
                   -z*L(L1n,:).*exp(-L(L1n,:)))./L(L1n,:).^2 ...
          + J(j2n,:).*(exp(-L(L1n,:)*z) - 1 ...
                   +z*L(L1n,:))./L(L1n,:).^2 ...
          + J(j3n,:).*(exp( L(L2n,:)*(z-1))-exp(-L(L2n,:)) ...
                   -z*L(L2n,:).*exp(-L(L2n,:)))./L(L2n,:).^2 ...
          + J(j4n,:).*(exp(-L(L2n,:)*z) - 1 ...
                   +z*L(L2n,:))./L(L2n,:).^2)/sigman;
      aux.names{indNeg(k)} = 'negPhis';
      aux.xLoc(indNeg(k)) = z;                                  
    end
  end

  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      phisTF(indPos(k),:) = (z ...
          +J(j1p,:).*(exp( L(L1p,:)*(z-1))-exp(-L(L1p,:)) ...
                  - z*L(L1p,:).*exp(-L(L1p,:)))./L(L1p,:).^2 ...
          +J(j2p,:).*(exp(-L(L1p,:)*z) - 1 ...
                  + z*L(L1p,:))./L(L1p,:).^2 ...
          +J(j3p,:).*(exp( L(L2p,:)*(z-1))-exp(-L(L2p,:)) ...
                  - z*L(L2p,:).*exp(-L(L2p,:)))./L(L2p,:).^2 ...
          +J(j4p,:).*(exp(-L(L2p,:)*z) - 1 ...
                  + z*L(L2p,:))./L(L2p,:).^2)/sigmap;
      aux.names{indPos(k)} = 'posPhis';
      aux.xLoc(indPos(k)) = 3-z;                                          
    end
  end

  if ~isempty(find(s==0,1))
      phisTF(:,s==0) = aux.dcGain;
  end

  % Computes DC gains (f->0)
  function dcGains = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize

    indNeg = find(locs <= 1);
    if ~isempty(indNeg)
      for kk = 1:length(indNeg)
        z = locs(indNeg(kk));
        dcGains(indNeg(kk)) = z*(z-2)/(2*sigman);
      end
    end

    indPos = find(locs >= 2);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        dcGains(indPos(kk)) = -z*(z-2)/(2*sigmap);
      end
    end
  end

  % Compute high-frequency gains (f->inf)
  function dTerm = calcDterm(locs)
    dTerm = 0*locs(:); % dummy initialize

    indNeg = find(locs <= 1);
    if ~isempty(indNeg)
      % set up negative-electrode constants...
      kappan = cellData.neg.kappa;
      Rfn = cellData.neg.Rf;
      Rdln = cellData.neg.Rdl;
      Rctn = Rct(1);
      ZseInfn = Rfn + 1./(1./Rctn + 1./Rdln);
      for kk = 1:length(indNeg)
        z = locs(indNeg(kk));
        if ZseInfn == 0
          dTerm(indNeg(kk)) = ((z-2)*z)/(2*sigman);
        else
          nuInfn = sqrt((1/sigman+1/kappan)/ZseInfn);
          dTerm(indNeg(kk)) = (-kappan*(cosh(nuInfn)-cosh(nuInfn*(z-1)))...
              -sigman*(1-cosh(nuInfn*z)+z*sinh(nuInfn)*nuInfn))...
              /(sigman*(sigman+kappan)*sinh(nuInfn)*nuInfn);
        end
      end
    end

    indPos = find(locs >= 2);
    if ~isempty(indPos)
      % set up positive-electrode constants...
      kappap = cellData.pos.kappa;
      Rfp = cellData.pos.Rf;
      Rdlp = cellData.pos.Rdl;
      Rctp = Rct(2);
      ZseInfp = Rfp + 1./(1./Rctp + 1./Rdlp);
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        if ZseInfp == 0
          dTerm(indPos(kk)) = -((3 - 4*(3-z) + (3-z)^2)/(2*sigmap));
        else
          nuInfp = sqrt((1/sigmap+1/kappap)/ZseInfp);
          dTerm(indPos(kk)) = (kappap*(cosh(nuInfp) - cosh(nuInfp*(z-1))) +...
              sigmap*(1-cosh(z*nuInfp)+z*nuInfp*sinh(nuInfp)))...
              /(sigmap*(kappap+sigmap)*nuInfp*sinh(nuInfp));
        end
      end
    end
  end
end