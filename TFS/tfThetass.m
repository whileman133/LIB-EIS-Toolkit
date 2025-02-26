% function [thetassTF,aux] = tfThetass(s,locs,cellData) 
% 
% Inputs:
%   s         = vector of 1j*w where w is vector of frequencies at which
%               the transfer function (TF) is to be evaluated
%   locs      = vector of locations in normalized coordinates (0..3) at
%               which the TF is to be evaluated
%   cellData  = data structure containing all parameter values for a cell
%               at a given setpoint, most likely the output of
%               evalSetpoint.m 
% Outputs:
%   thetassTF = a matrix of the TF evaluated at every combination of "s"
%               and "locs" -- each row holds the TF for one location and
%               all frequencies; there is one row for every location
%   aux       = other outputs sometimes needed by calling routines (e.g.,
%               the xRAs). These outputs include dc-gains, high-frequency
%               gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the INTEGRATOR-REMOVED TF for normalized
% concentration of lithium in the solid, at the surface. If you wish to
% evaluate the TF with the integrator, use tfThetassInt.m instead.
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

function [thetassTF,aux] = tfThetass(s,locs,cellData) 
  % Initialize thetassstar to thetass
  s = s(:).';          % force "s" to be a row vector
  [C,Lambda,J,Z,Rct] = tfCommon(s,cellData); % Get Rct; store for efficiency
  cellData.common.C = C; cellData.common.L = Lambda; 
  cellData.common.J = J; cellData.common.Z = Z; cellData.common.Rct = Rct;
  [thetassTF,~] = tfThetassInt(s,locs,cellData);

  % Initialize outputs
  [aux.dcGain,r0n,r0p] = dcGain(locs);
  aux.hfGain = 0*locs(:);                % high-frequency gain  
  aux.res0 = 0*locs(:);                  % integrator residue
  aux.names = cell(length(locs),1);      % TF names
  aux.xLoc = zeros(1,length(locs));      % TF locations  


  indNeg = find(locs <= 1);
  if ~isempty(indNeg)
    aux.res0(indNeg) = r0n;
    thetassTF(indNeg,:) = thetassTF(indNeg,:) ...
                              - r0n*ones(size(indNeg(:)))*(1./s);
    [aux.names{indNeg}] = deal('negThetass');
    aux.xLoc(indNeg) = locs(indNeg);                                                                        
  end  

  indPos = find(locs >= 2);
  if ~isempty(indPos)
    aux.res0(indPos) = r0p;
    thetassTF(indPos,:) = thetassTF(indPos,:) ...
                              - r0p*ones(size(indPos(:)))*(1./s);
    [aux.names{indPos}] = deal('posThetass');
    aux.xLoc(indPos) = locs(indPos);                                                                                                    
  end
  inds0 = find(s==0,1);
  if ~isempty(inds0)
    thetassTF(:,inds0) = aux.dcGain;
  end

  function [dcGains,r0n,r0p] = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize
    
    T = cellData.const.T;
    
    % Set up constant values
    psi = cellData.const.psi;
    kappaD = cellData.const.kD;
    Q = cellData.const.Q; 
    
    % set up negative-electrode constants...
    Dsn = cellData.neg.Ds;
    DeltaQn = abs(cellData.neg.theta100 - cellData.neg.theta0);
    csmaxn = 10800*Q*Dsn/DeltaQn;
    csmaxn2 = 3*Q*Dsn/DeltaQn;
    dUdqn = cellData.neg.dUocp;
    Cdln = cellData.neg.Cdl;
    sigman = cellData.neg.sigma;
    kappan = cellData.neg.kappa;
    qen = cellData.neg.qe;
    Rctn = cellData.common.Rct(1);
    nDLn = cellData.neg.nDL;    
    Rdln = cellData.neg.Rdl;
    wDLn = cellData.neg.wDL;
    CdlEffn = Cdln^(2-nDLn)*wDLn^(nDLn-1);
    r0n = 3*Dsn/(3*CdlEffn*Dsn*dUdqn-csmaxn);

    % set up separator constants
    kappas = cellData.sep.kappa; % not SOC dependent
    qes = cellData.sep.qe;       % not SOC dependent

    % set up positive-electrode constants...
    Dsp = cellData.pos.Ds;
    DeltaQp = abs(cellData.pos.theta100 - cellData.pos.theta0);
    csmaxp = 10800*Q*Dsp/DeltaQp;
    csmaxp2 = 3*Q*Dsp/DeltaQp;
    dUdqp = cellData.pos.dUocp;
    Cdlp = cellData.pos.Cdl;
    sigmap = cellData.pos.sigma;
    kappap = cellData.pos.kappa;
    qep = cellData.pos.qe;
    Rctp = cellData.common.Rct(2);
    nDLp = cellData.pos.nDL;    
    Rdlp = cellData.pos.Rdl;
    wDLp = cellData.pos.wDL;
    CdlEffp = Cdlp^(2-nDLp)*wDLp^(nDLp-1);    
    r0p = -3*Dsp/(3*CdlEffp*Dsp*dUdqp-csmaxp); 

    % Zsestar as s-> 0 
    if CdlEffn == 0
      z0n = -1/(5*csmaxn);
    else
      z0n = (csmaxn*(-1+15*CdlEffn*Dsn*Rctn) - ...
             45*CdlEffn*Dsn^2*dUdqn*((nDLn-1)*Cdln/wDLn - ...
             CdlEffn*Rdln))/(5*(-3*CdlEffn*Dsn*dUdqn + csmaxn)^2);
    end    
    Qpp1den = 18*Dsn*dUdqn*(qen+qep+qes)*kappap* ...
              kappan^2*kappas*sigman*(psi*T)^2;            
    indNeg = find(locs <= 1);
    if ~isempty(indNeg)
      for kkk = 1:length(indNeg)
        z = locs(indNeg(kkk));
        Qe0 = (qen*kappap*kappas... % dc gain of Thetae(z,0)
            +qep*(2*kappan*kappas+3*kappap*kappas+6*kappan*kappap)...
            +qes*(3*kappap*kappas+3*kappan*kappap))...
            /(6*kappan*kappap*kappas*(psi*T)*(qen+qep+qes))...
            -z^2/(2*kappan*(psi*T)); 
        Qpp1num = -(-3*Dsn*dUdqn*qen*(qen*(1-3*z^2)*kappap*kappas + ...
                   2*qep*kappan*kappas + ...
                   3*qes*kappap*(kappan + kappas - z^2*kappas) + ...
                   3*qep*kappap*(2*kappan + kappas - ...
                   z^2*kappas))*sigman + ...
                   csmaxn2*(qen + qep + qes)*kappap*kappas*(T*(-1+3*z^2)*kappaD*sigman + ...
                   ((-2+6*z-3*z^2)*kappan+sigman - ...
                   3*z^2*sigman)*(psi*T)));
          
        dcGains(indNeg(kkk)) = z0n + r0n*3600*(qen*Qe0 - (psi*T)*kappan*Qpp1num/Qpp1den);
      end
    end

    % Zsestar as s-> 0 
    if CdlEffp == 0
      z0p = -1/(5*csmaxp);
    else
      z0p = (csmaxp*(-1+15*CdlEffp*Dsp*Rctp) - ...
             45*CdlEffp*Dsp^2*dUdqp*((nDLp-1)*Cdlp/wDLp - ...
             CdlEffp*Rdlp))/(5*(-3*CdlEffp*Dsp*dUdqp + csmaxp)^2);
    end    
    Qpp1den = 18*Dsp*dUdqp*(qen+qep+qes)*kappan* ...
              kappap^2*kappas*sigmap*(psi*T)^2;
    indPos = find(locs >=2);
    if ~isempty(indPos)
      for kkk = 1:length(indPos)
        z = 3-locs(indPos(kkk));
        Qe0 = (qen*(-2*kappap*kappas-3*kappan*kappas-6*kappan*kappap)...
            -qep*kappan*kappas+qes*(-3*kappan*kappas-3*kappan*kappap))...
            /(6*kappan*kappap*kappas*(psi*T)*(qen+qep+qes))...
            +z^2/(2*kappap*(psi*T));
        Qpp1num = (-3*Dsp*dUdqp*qep*(qep*(1-3*z^2)*kappan*kappas + ...
                   2*qen*kappap*kappas + ...
                   3*qes*kappan*(kappap + kappas - z^2*kappas) + ...
                   3*qen*kappan*(2*kappap + kappas - ...
                   z^2*kappas))*sigmap + ...
                   csmaxp2*(qen + qep + qes)*kappan*kappas*(T*(-1+3*z^2)*kappaD*sigmap + ...
                   ((-2+6*z-3*z^2)*kappap+sigmap - ...
                   3*z^2*sigmap)*(psi*T)));
        dcGains(indPos(kkk)) = -(z0p + r0p*3600*(qep*Qe0 - (psi*T)*kappap*Qpp1num/Qpp1den));
      end
    end    
  end
end  