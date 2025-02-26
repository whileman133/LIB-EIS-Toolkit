% function [phieTF,aux] = tfPhie(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   phieTF   = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the electrolyte-potential TF.
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

function [phieTF,aux] = tfPhie(s,locs,cellData)
  % Set up some constants
  T = cellData.const.T;
    
  % Initialize shared variables needed for calculating gains
  kappan = cellData.neg.kappa;
  kappap = cellData.pos.kappa;
  kappas = cellData.sep.kappa;
  kappaD = cellData.const.kD;
  psi = cellData.const.psi;

  % Initialize outputs
  s = s(:).';          % force "s" to be a row vector
  [C,L,J,Z,Rct] = tfCommon(s,cellData); 
  cellData.common.C = C; cellData.common.L = L; 
  cellData.common.J = J; cellData.common.Z = Z;
  cellData.common.Rct = Rct;

  phieTF = zeros(length(locs),length(s));
  aux.hfGain = calcDterm(locs);          % high-frequency gain
  aux.dcGain = dcGain(locs);             % low-frequency gain
  aux.res0 = 0*locs(:);                  % integrator residue
  aux.names = cell(length(locs),1);      % TF names
  aux.xLoc = zeros(1,length(locs));      % TF locations

  % Get common components... (next three lines needed to keep MATLAB happy)
  L1n=[];L2n=[];L1s=[];L1p=[];L2p=[];c1n=[];c2n=[];c3n=[];c4n=[];c1s=[];
  c2s=[];c1p=[];c2p=[];c3p=[];c4p=[];j1n=[];j2n=[];j3n=[];j4n=[];
  j1p=[];j2p=[];j3p=[];j4p=[];Zsen=[];Zsep=[];Zsn=[];Zsp=[];Isn=[];Isp=[]; %#ok<NASGU>
  eval(cellData.common.ind); % set c##,L##,J##,Z## indices

  % Compute TF at each electrode
  % [phien]1 & [phisen]2
  indNeg = find(locs <= 1);
  if ~isempty(indNeg)
    for k = 1:length(indNeg)
      z = locs(indNeg(k));
      phieTF(indNeg(k),:) = -( ...
          +J(j1n,:).*(exp( L(L1n,:)*(z-1))-exp(-L(L1n,:)) ...
                   - z*L(L1n,:).*exp(-L(L1n,:)))./L(L1n,:).^2 ...
          +J(j2n,:).*(exp(-L(L1n,:)*z) - 1 ...
                   + z*L(L1n,:))./L(L1n,:).^2 ...
          +J(j3n,:).*(exp( L(L2n,:)*(z-1))-exp(-L(L2n,:)) ...
                   - z*L(L2n,:).*exp(-L(L2n,:)))./L(L2n,:).^2 ...
          +J(j4n,:).*(exp(-L(L2n,:)*z) - 1 ...
                   + z*L(L2n,:))./L(L2n,:).^2)/kappan;
      phieTF(indNeg(k),:) = phieTF(indNeg(k),:) ...
          -kappaD*T*(C(c1n,:).*(exp(L(L1n,:)*(z-1))-exp(-L(L1n,:))) ...
                   + C(c2n,:).*(exp(-L(L1n,:)*z)-1) ...
                   + C(c3n,:).*(exp(L(L2n,:)*(z-1))-exp(-L(L2n,:))) ...
                   + C(c4n,:).*(exp(-L(L2n,:)*z)-1));
      aux.names{indNeg(k)} = 'negPhie';
      aux.xLoc(indNeg(k)) = z;                          
    end
  end
  z = 1;
  phien1 = -( ...
      +J(j1n,:).*(exp( L(L1n,:)*(z-1))-exp(-L(L1n,:)) ...
               - z*L(L1n,:).*exp(-L(L1n,:)))./L(L1n,:).^2 ...
      +J(j2n,:).*(exp(-L(L1n,:)*z) - 1 ...
               + z*L(L1n,:))./L(L1n,:).^2 ...
      +J(j3n,:).*(exp( L(L2n,:)*(z-1))-exp(-L(L2n,:)) ...
               - z*L(L2n,:).*exp(-L(L2n,:)))./L(L2n,:).^2 ...
      +J(j4n,:).*(exp(-L(L2n,:)*z) - 1 ...
               + z*L(L2n,:))./L(L2n,:).^2)/kappan;
  phien1 = phien1 ...
      -kappaD*T*(C(c1n,:).*(exp(L(L1n,:)*(z-1))-exp(-L(L1n,:))) ...
               + C(c2n,:).*(exp(-L(L1n,:)*z)-1) ...
               + C(c3n,:).*(exp(L(L2n,:)*(z-1))-exp(-L(L2n,:))) ...
               + C(c4n,:).*(exp(-L(L2n,:)*z)-1));

  % [phies]1 & [phies]2
  indSep = find(locs > 1 & locs < 2);
  if ~isempty(indSep)
    for k = 1:length(indSep)
      z = locs(indSep(k))-1;
      phieTF(indSep(k),:) = phien1 - z/kappas;
      phieTF(indSep(k),:) = phieTF(indSep(k),:) ...
          -kappaD*T*(C(c1s,:).*exp(L(L1s,:)*(z-1)) - C(c1s,:).*exp(-L(L1s,:)) ...
                   + C(c2s,:).*(exp(-L(L1s,:)*z)-1));
      aux.names{indSep(k)} = 'sepPhie';
      aux.xLoc(indSep(k)) = z+1;        
    end
  end
  z = 1;
  phies1 = phien1 - z/kappas;
  phies1 = phies1 ...
      -kappaD*T*(C(c1s,:).*exp(L(L1s,:)*(z-1)) - C(c1s,:).*exp(-L(L1s,:)) ...
               + C(c2s,:).*(exp(-L(L1s,:)*z)-1));

  % [phiep]1 & [phiep]2
  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      phieTF(indPos(k),:) = phies1  + (z-1)/kappap - ( ...
          +J(j1p,:).*(exp( L(L1p,:)*(z-1)) - 1           +(1-z)*L(L1p,:)                )./L(L1p,:).^2 ...
          +J(j2p,:).*(exp(-L(L1p,:)*z)    -exp(-L(L1p,:))-(1-z)*L(L1p,:).*exp(-L(L1p,:)))./L(L1p,:).^2 ...
          +J(j3p,:).*(exp( L(L2p,:)*(z-1))-1             +(1-z)*L(L2p,:)                )./L(L2p,:).^2 ...
          +J(j4p,:).*(exp(-L(L2p,:)*z)    -exp(-L(L2p,:))-(1-z)*L(L2p,:).*exp(-L(L2p,:)))./L(L2p,:).^2)/kappap;
      phieTF(indPos(k),:) = phieTF(indPos(k),:) ...
          -kappaD*T*(C(c1p,:).*(exp( L(L1p,:)*(z-1))-1             ) + ...
                     C(c2p,:).*(exp(-L(L1p,:)*z)    -exp(-L(L1p,:))) + ...
                     C(c3p,:).*(exp( L(L2p,:)*(z-1))-1             ) + ...
                     C(c4p,:).*(exp(-L(L2p,:)*z)    -exp(-L(L2p,:))));
      aux.names{indPos(k)} = 'posPhie';
      aux.xLoc(indPos(k)) = 3-z;                           
    end
  end
  if ~isempty(find(s==0,1))
      phieTF(:,s==0) = dcGain(locs);
  end

  % Computes DC gains (f->0)
  function dcGains = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize

    indNeg = find(locs <= 1);
    if ~isempty(indNeg)
      for kk = 1:length(indNeg)
        z = locs(indNeg(kk));
        dcGains(indNeg(kk)) = -z^2/(2*kappan)+kappaD*T*z^2/(2*(psi*T)*kappan);
      end
    end

    dcNeg1 = -1/(2*kappan) + kappaD*T/(2*(psi*T)*kappan);
    indSep = find(locs > 1 & locs < 2);
    if ~isempty(indSep)
      for kk = 1:length(indSep)
        z = locs(indSep(kk))-1;
        dcGains(indSep(kk)) = dcNeg1-z/kappas+kappaD*T*z/((psi*T)*kappas);
      end
    end

    dcSep1 = dcNeg1 - 1/kappas + kappaD*T/((psi*T)*kappas);
    indPos = find(locs >= 2);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        dcGains(indPos(kk)) = dcSep1 + (z^2-1)/(2*kappap) ...
            - (z^2-1)*kappaD*T/(2*(psi*T)*kappap);
      end
    end
  end

 % A function computes high-frequency gains (f->inf)
  function dTerm = calcDterm(locs)
    dTerm = 0*locs(:); % dummy initialize

    Rf = cellData.neg.Rf;
    Rctn = cellData.common.Rct(1);
    ZseInfn = Rf + 1./(1./Rctn + 1./cellData.neg.Rdl);
    sigman = cellData.neg.sigma;

    Rf = cellData.pos.Rf;
    Rctp = cellData.common.Rct(2);
    ZseInfp = Rf + 1./(1./Rctp + 1./cellData.pos.Rdl);
    sigmap = cellData.pos.sigma;

    indNeg = find(locs <= 1);
    nuInfn = sqrt((1/sigman+1/kappan)/ZseInfn);
    if ~isempty(indNeg)
      for kk = 1:length(indNeg)
        z = locs(indNeg(kk));
        if ZseInfn == 0
          dTerm(indNeg(kk)) = -(z^2/(2*kappan));
        else
          dTerm(indNeg(kk)) = (sigman*(1-cosh(nuInfn*z))...
              +kappan*(cosh(nuInfn))...
              -kappan*cosh(nuInfn*(1-z))...
              -kappan*z*sinh(nuInfn)*nuInfn)...
              /(kappan*(sigman+kappan)*sinh(nuInfn)*nuInfn);
        end
      end
    end

    indSep = find(locs > 1 & locs < 2);
    if ~isempty(indSep)
      for kk = 1:length(indSep)
        z = 1;
        if nuInfn == 0
          dTerm(indSep(kk)) = -1/(2*kappan);
        else
          dTerm(indSep(kk)) = (sigman*(1-cosh(nuInfn*z))...
              +kappan*(cosh(nuInfn))...
              -kappan*cosh(nuInfn*(1-z))...
              -kappan*z*sinh(nuInfn)*nuInfn)...
              /(kappan*(sigman+kappan)*sinh(nuInfn)*nuInfn);
        end

        z = locs(indSep(kk))-1;
        dTerm(indSep(kk)) = dTerm(indSep(kk)) -z/kappas;
      end
    end

    indPos = find(locs >= 2);
    if ~isempty(indPos)
      nuInfp = sqrt((1/sigmap+1/kappap)/ZseInfp);
      for kk = 1:length(indPos)
        z = 1;
        if nuInfn == 0
          dTerm(indPos(kk)) = -1/(2*kappan) - 1/kappas;
        else
          dTerm(indPos(kk)) = (sigman*(1-cosh(nuInfn*z))...
              +kappan*(cosh(nuInfn))...
              -kappan*cosh(nuInfn*(1-z))...
              -kappan*z*sinh(nuInfn)*nuInfn)...
              /(kappan*(sigman+kappan)*sinh(nuInfn)*nuInfn) - 1/kappas;
        end

        z = 3-locs(indPos(kk));
        if nuInfp == 0
          dTerm(indPos(kk)) = dTerm(indPos(kk)) ...
            -(((1-z)*(sigmap*(-z-1)-kappap*(3-z)))/...
            (2*kappap*(kappap+sigmap)));
        else
          dTerm(indPos(kk)) = dTerm(indPos(kk)) ...
              -(sigmap*(cosh(nuInfp)-cosh(nuInfp*z))...
              +kappap*(1-cosh(nuInfp*(1-z))+(1-z)*sinh(nuInfp)*nuInfp))...
              /(kappap*(sigmap+kappap)*sinh(nuInfp)*nuInfp);
        end
      end
    end
  end
end