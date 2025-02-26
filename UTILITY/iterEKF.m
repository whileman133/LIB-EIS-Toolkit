% function [zk,boundzk,ekfData,Xind] = iterEKF(vk,ik,Tk,ekfData)
% 
% Inputs:
%   vk       = The measured cell voltage for this time sample (V)
%   ik       = The measured cell current for this time sample (A)
%   Tk       = The measured cell temperature for this time sample (degC)
%   ekfData  = Data structure initialized with initKF and updated here
% Outputs:
%   zk       = Comprises three sections, stacked vertically. The top of zk
%              holds all output variables (with nonlinear corrections
%              applied) in the same order they are organized in the ROMs.
%              Under this, zk has a (scalar) voltage estimate. Under this,
%              zk has a (scalar) SOC estimate (0..1).
%   boundzk  = A vector of 3-sigma confidence bounds for every element in
%              zk. Same dimension as zk.
%   ekfData  = Data structure used by EKF, with updated contents.
%   Xind     = Optional output, containing indices and scaling factors for
%              the individual ROMs used in the blend (this should not be
%              needed by most users so is a utility/debug output only).
%
% This function updates the EKF based on measurements made during one
% sampling period. It outputs estimates of internal variables, voltage, and
% SOC, as well as 3-sigma confidence bounds on those estimates.
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

function [zk,boundzk,ekfData,Xind] = iterEKF(vk,ik,Tk,ekfData)
  persistent ind loc numT numZ Tpts Zpts n warnCount cellData F R
  if isempty(ind)
    setupIndsLocs;
    % All temperature and SOC setpoints (sort in ascending order)
    Tpts = unique(ekfData.T);
    Zpts = unique(ekfData.Z);
    numT = length(Tpts); numZ = length(Zpts);
    n = ekfData.n;
    ekfData.ind = ind; % helpful to save for later when plotting
    ekfData.loc = loc;
    ekfData.Q = ekfData.cellData.function.const.Q(); % cell capacity in Ah    
    warnCount = 0;
    cellData  = ekfData.cellData;
    F         = ekfData.cellData.const.F;
    R         = ekfData.cellData.const.R;     
    
    if isempty(ekfData.SOC0)
      ocvFn = @(x) vk - (cellData.function.pos.Uocp(cellData.function.pos.soc(x)) ...
                 - cellData.function.neg.Uocp(cellData.function.neg.soc(x)));
      ekfData.SOC0 = fzero(ocvFn,0.5); % guess at initial SOC by voltage 
      fprintf('Initializing SOC0 to %g%%\n',ekfData.SOC0*100);
    end
  end

  if warnCount > ekfData.maxWarn % EKF is probably broken/lost
    zk = NaN(ekfData.nz+2,1); boundzk = zk;    
    Xind.gamma = NaN(4,1); Xind.theT = NaN(4,1); Xind.theZ = NaN(4,1);
    return
  end
  
  % Convert degC to K
  if Tk > 100
    warning('iterSPKF assumes that Tk is in Celsius. Converting...');
  else
    Tk = Tk + 273.15; % Convert to Kelvin...
  end
  
  % ----------------------------------------------------------------------
  % Steps 1a and 1b... update state(s) prediction and coviariance
  switch ekfData.method
    case 'OB'
      % Need to update state and covariance for ALL models
      for theT = 1:numT
        for theZ = 1:numZ
          A = ekfData.M(theT,theZ).A;
          xhat = ekfData.M(theT,theZ).xhat;
          xhat = A.*xhat + ekfData.priorI; % don't forget to use prior current!
          ekfData.M(theT,theZ).xhat = xhat;

          SigmaX = ekfData.M(theT,theZ).SigmaX;
          SigmaX = diag(A)*SigmaX*diag(A) + ekfData.SigmaW;
          ekfData.M(theT,theZ).SigmaX = SigmaX;
        end
      end
      ekfData.x0      = ekfData.x0 + ekfData.priorI; % not ik!!
      ekfData.SigmaX0 = ekfData.SigmaX0 + ekfData.SigmaW;

      % Compute time-update prediction of SOC
      SOC = ekfData.SOC0 - ekfData.x0*(ekfData.Ts/(3600*ekfData.Q)); 
    case 'MB'
      % Need to update state and covariance of single model
      SOC = ekfData.SOC0 - ekfData.xhat(end)*(ekfData.Ts/(3600*ekfData.Q)); 
      Xind = getXind(Tk,SOC);
      A1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).A;
      A2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).A;
      A3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).A;
      A4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).A;
      AMB = [[A1,A2,A3,A4]*Xind.gamma;1]; % model-blend A, with integrator
      ekfData.xhat = AMB.*ekfData.xhat + ekfData.priorI;
      
      ekfData.SigmaX = diag(AMB)*ekfData.SigmaX*diag(AMB) + ekfData.SigmaW;
      SOC = ekfData.SOC0 - ekfData.xhat(end)*(ekfData.Ts/(3600*ekfData.Q)); 
  end

  % ----------------------------------------------------------------------
  % Step 1c... predict measured voltage
  Xind = getXind(Tk,SOC);
  [vhat,Z] = getVariables(ik,Xind,Tk); % this uses ik, not priorI
  
  % ----------------------------------------------------------------------
  % Step 2a... operates only on four nearest neighbors
  switch ekfData.method
    case 'OB'
      [ChatV,C0] = getChatV(Xind,Z,Tk);
      SigmaVtilde = cell(4,1); L = SigmaVtilde; % reserve space
      for theModel = 1:4
        SigmaX = ekfData.M(Xind.theT(1),Xind.theZ(1)).SigmaX;
        SigmaVtilde{theModel} = ChatV{theModel}*SigmaX*ChatV{theModel}' ...
                                 + ekfData.SigmaV;
        L{theModel} = SigmaX*ChatV{theModel}'/SigmaVtilde{theModel};
      end
      SigmaX0 = ekfData.SigmaX0;
      SigmaVtilde0 = (C0*SigmaX0*C0' + ekfData.SigmaV);
      L0 = SigmaX0*C0/SigmaVtilde0;
    case 'MB'
      ChatV = getChatV(Xind,Z,Tk);
      SigmaVtilde = ChatV*ekfData.SigmaX*ChatV' + ekfData.SigmaV;
      L = ekfData.SigmaX*ChatV'/SigmaVtilde;
  end
  
  % ----------------------------------------------------------------------
  % Steps 2b and 2c... operate only on four nearest neighbors
  residual = vk - vhat;
  switch ekfData.method
    case 'OB'
      for theModel = 1:4
        xhat = ekfData.M(Xind.theT(theModel),Xind.theZ(theModel)).xhat;
        xhat = xhat + L{theModel}*residual;
        ekfData.M(Xind.theT(theModel),Xind.theZ(theModel)).xhat = xhat;

        SigmaX = ekfData.M(Xind.theT(theModel),Xind.theZ(theModel)).SigmaX;
        SigmaX = SigmaX - L{theModel}*SigmaVtilde{theModel}*(L{theModel})';
        [~,SS,VV] = svd(SigmaX);
        HH = VV*SS*VV';
        SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness

        % Q-bump code
        if residual^2>9*SigmaVtilde{theModel} % bad voltage estimate by 3-SigmaX, bump Q 
          fprintf('Increasing SigmaX\n');
          SigmaX = SigmaX*2;
        end

        ekfData.M(Xind.theT(theModel),Xind.theZ(theModel)).SigmaX = SigmaX;
      end
      ekfData.x0 = ekfData.x0 + L0*residual;
      ekfData.SigmaX0 = ekfData.SigmaX0 - L0*SigmaVtilde0*L0;
      
      % Compute measurement-update estimate of SOC
      SOC = ekfData.SOC0 - ekfData.x0*(ekfData.Ts/(3600*ekfData.Q));         
    case 'MB'
      ekfData.xhat = ekfData.xhat + L*residual;
      SigmaX = ekfData.SigmaX;
      SigmaX = SigmaX - L*SigmaVtilde*L';
      [~,SS,VV] = svd(SigmaX);
      HH = VV*SS*VV';
      SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness

      % Q-bump code
      if residual^2>9*SigmaVtilde % bad voltage estimate by 3-SigmaX, bump Q 
        fprintf('Increasing SigmaX\n');
        SigmaX = SigmaX*2;
      end
      ekfData.SigmaX = SigmaX;      

      % Compute measurement-update estimate of SOC
      SOC = ekfData.SOC0 - ekfData.xhat(end)*(ekfData.Ts/(3600*ekfData.Q));         
  end
      
  % ----------------------------------------------------------------------
  % Steps 3a and 3b... operate only on four nearest neighbors
  Xind = getXind(Tk,SOC);
  [vhat,zk,Zsoc] = getVariables(ik,Xind,Tk); 
  zk = [zk;vhat;Zsoc];
    
  res = -ekfData.Ts/(3600*ekfData.Q);
  switch ekfData.method
    case 'OB'
      SigmaZ = zeros(ekfData.nz); % initialize 
      SigmaV = 0;
      [ChatZ,ChatV,ChatZ0,ChatV0] = getChatZ(Xind,zk,Tk);
      for theModel = 1:4
        SigmaX = ekfData.M(Xind.theT(1),Xind.theZ(1)).SigmaX;
        SigmaZ = SigmaZ + ChatZ{theModel}*SigmaX*ChatZ{theModel}';
        SigmaV = SigmaV + ChatV{theModel}*SigmaX*ChatV{theModel}';
      end
      SigmaZ = SigmaZ + ChatZ0*ekfData.SigmaX0*ChatZ0';
      SigmaV = SigmaV + ChatV0*ekfData.SigmaX0*ChatV0';
      SigmaSOC = res*ekfData.SigmaX0*res;
    case 'MB'
      [ChatZ,ChatV] = getChatZ(Xind,zk,Tk);
      SigmaZ = ChatZ*ekfData.SigmaX*ChatZ';
      SigmaV = ChatV*ekfData.SigmaX*ChatV';
      SigmaSOC = res*ekfData.SigmaX(end,end)*res;
  end
  boundzk = 3*sqrt([diag(SigmaZ); SigmaV; SigmaSOC]);
  
  
  % ----------------------------------------------------------------------
  % Time to return 
  ekfData.priorI = ik;

  %% ======================================================================
  % The functions below this point implement the details of the higher-
  % level functionality indicated above
  % =======================================================================
  
  %% ----------------------------------------------------------------------
  % This function computes the indexing variables for the four blended mdls
  function Xind = getXind(Tk,SOC)
    % Find the two closest Zspts setpoints: "Zupper" and "Zlower"
    dZ = abs(SOC-Zpts); [~,iZ] = sort(dZ);
    Zupper = Zpts; iZupper = iZ; % default for single-model ROM
    Zlower = Zpts; iZlower = iZ; % default for single-model ROM
    if length(iZ)>1
      Zupper = Zpts(iZ(1)); iZupper = iZ(1);
      Zlower = Zpts(iZ(2)); iZlower = iZ(2);
      if Zupper<Zlower
        Zupper = Zpts(iZ(2)); iZupper = iZ(2);
        Zlower = Zpts(iZ(1)); iZlower = iZ(1);
      end
    end
    
    % Find the two closest temperature setpoints: "Tupper" and "Tlower"
    dT = abs(Tk-Tpts); [~,iT] = sort(dT);
    Tupper = Tpts; iTupper = iT; % default for single-model ROM
    Tlower = Tpts; iTlower = iT; % default for single-model ROM
    if length(iT)>1
      Tupper = Tpts(iT(1)); iTupper = iT(1);
      Tlower = Tpts(iT(2)); iTlower = iT(2);
      if Tupper<Tlower
        Tupper = Tpts(iT(2)); iTupper = iT(2);
        Tlower = Tpts(iT(1)); iTlower = iT(1);
      end
    end
    
    alphaZ = 0; alphaT = 0;
    if length(iZ)>1, alphaZ = (SOC-Zlower)/(Zupper-Zlower); end
    if length(iT)>1, alphaT = (Tk-Tlower)/(Tupper-Tlower); end    
    Xind.gamma = [(1-alphaT)*(1-alphaZ); % 00
                  (1-alphaT)*alphaZ;     % 01
                  alphaT*(1-alphaZ);     % 10
                  alphaT*alphaZ];        % 11
    Xind.theT = [iTlower iTlower iTupper iTupper]';
    Xind.theZ = [iZlower iZupper iZlower iZupper]';
  end

  %% ----------------------------------------------------------------------
  % This function computes all cell nonlinear output variables
  function [Vcell,Z,Zsoc] = getVariables(ik,Xind,T)

    % ---------------------------------------------------------------------
    % Step 1: Extract the state sigma points of 4 nearest-neighbor models
    %         (needed for OB only)
    % ---------------------------------------------------------------------
    switch ekfData.method
      case 'OB'
        xk1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).xhat;
        xk2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).xhat;
        xk3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).xhat;
        xk4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).xhat;
        x0  = ekfData.x0;
      case 'MB'
        xhat = ekfData.xhat;
        x0 = xhat(end);
    end
    % Compute SOC
    xSOC = ekfData.SOC0 - x0*(ekfData.Ts/(3600*ekfData.Q)); 

    % ---------------------------------------------------------------------
    % Step 2: Initialize some variables/constants we will need
    % ---------------------------------------------------------------------
    SOCnAvg = cellData.function.neg.soc(xSOC,T);
    SOCpAvg = cellData.function.pos.soc(xSOC,T);
    if SOCnAvg < 0
      shortWarn('SOCnAvg < 0'); SOCnAvg = 1e-6;
    end
    if SOCnAvg > 1
      shortWarn('SOCnAvg > 1'); SOCnAvg = 1-1e-6;
    end
    if SOCpAvg < 0
      shortWarn('SOCpAvg < 0'); SOCpAvg = 1e-6;
    end
    if SOCpAvg > 0.998
      shortWarn('SOCpAvg > 1'); SOCpAvg = 0.998;
    end
    
    % ---------------------------------------------------------------------
    % Step 3: Find the linear outputs as
    %         Z[k] = C*x[k] + D*iapp[k]
    % ---------------------------------------------------------------------
    Z = zeros(ekfData.nz,1);
    C1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).C;
    C2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).C;
    C3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).C;
    C4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).C;
    D1term = ekfData.M(Xind.theT(1),Xind.theZ(1)).D*ik;
    D2term = ekfData.M(Xind.theT(2),Xind.theZ(2)).D*ik;
    D3term = ekfData.M(Xind.theT(3),Xind.theZ(3)).D*ik;
    D4term = ekfData.M(Xind.theT(4),Xind.theZ(4)).D*ik;
    switch ekfData.method
      case 'OB'
        Z = [C1*xk1 + D1term, C2*xk2 + D2term, ...
                  C3*xk3 + D3term, C4*xk4 + D4term]*Xind.gamma;
      case 'MB'
        Z = [C1*xhat(1:n) + D1term, C2*xhat(1:n) + D2term, ...
                  C3*xhat(1:n) + D3term, C4*xhat(1:n) + D4term]*Xind.gamma;
    end
    
    % --------------------------------------------------------
    % Step 4: Apply nonlinear corrections to the linear output
    % --------------------------------------------------------
    % No corrections needed for ifdl, if, or idl
    % But, we need If at current collectors later on
    If0 = Z(ind.If0);
    If3 = Z(ind.If3);

    % Solid surface stoichiometries (thetass)
    % (we are guaranteed that ind.negThetass and ind.posThetass not empty)
    % Note that we are implementing integrator-removed Thetass, so we need
    % to add in SOCnAvg instead of SOC0n (and SOCpAvg instead of SOC0p)
    Z(ind.negThetass) = Z(ind.negThetass) + SOCnAvg;
    if any(Z(ind.negThetass) < 0)
      shortWarn('negThetass < 0'); 
      ZnegThetass = Z(ind.negThetass);
      ZnegThetass(ZnegThetass<0) = 1e-6;
      Z(ind.negThetass) = ZnegThetass;
    end
    if any(Z(ind.negThetass) > 1)
      shortWarn('negThetass > 1'); 
      ZnegThetass = Z(ind.negThetass);
      ZnegThetass(ZnegThetass>1) = 1-1e-6;
      Z(ind.negThetass) = ZnegThetass;
    end
    
    Z(ind.posThetass) = Z(ind.posThetass) + SOCpAvg;
    if any(Z(ind.posThetass) < 0)
      shortWarn('posThetass < 0'); 
      ZposThetass = Z(ind.posThetass);
      ZposThetass(ZposThetass<0) = 1e-6;
      Z(ind.posThetass) = ZposThetass;
    end
    if any(Z(ind.posThetass) > 0.998)
      shortWarn('posThetass > 1'); 
      ZposThetass = Z(ind.posThetass);
      ZposThetass(ZposThetass>0.998) = 0.998;
      Z(ind.posThetass) = ZposThetass;
    end

    % Solid-electrolyte potential difference (phise)
    % The linear output from Z is integrator-removed version 
    % (we are guaranteed that ind.negPhise not empty)
    UocpnAvg = cellData.function.neg.Uocp(SOCnAvg,T);
    UocppAvg = cellData.function.pos.Uocp(SOCpAvg,T);
    Z(ind.negPhise) = Z(ind.negPhise) + UocpnAvg;
    if ~isempty(ind.posPhise)
      Z(ind.posPhise) = Z(ind.posPhise) + UocppAvg;
    end

    % Compute electrolyte potential: first phie(0,t) then phie(1:3,t)
    % (we are guaranteed that ind.Phie not empty)
    PhieTilde3 = Z(ind.Phie(end));
    Phise0 = Z(ind.Phise0); 
    for theLoc = 1:length(ind.Phie)
      if loc.Phie(theLoc) == 0
        Z(ind.Phie(theLoc)) = 0 - Phise0;
      else
        Z(ind.Phie(theLoc)) = Z(ind.Phie(theLoc),:) - Phise0;
      end
    end

    % Compute electrolyte stoichiometries (thetae)
    % (we are guaranteed that ind.Thetae not empty)
    Z(ind.Thetae) = Z(ind.Thetae) + 1;
    if any(Z(ind.Thetae,:) < 0)
      shortWarn('Thetae < 0'); 
      ZThetae = Z(ind.Thetae);
      ZThetae(Zthetae < 0) = 1e-6;
      Z(ind,Thetae) = ZThetae;
    end

    % Compute overpotential at current-collectors via asinh method (eta)
    k0n = cellData.function.neg.k0(SOCnAvg,T);
    k0p = cellData.function.pos.k0(SOCpAvg,T);

    i0n = k0n*sqrt(Z(ind.Thetae(1)).* ...
               (1-Z(ind.Thetass0)).*Z(ind.Thetass0));
    i0p = k0p*sqrt(Z(ind.Thetae(end)).*...
               (1-Z(ind.Thetass3)).*Z(ind.Thetass3));

    negEta0 = 2*R*T/F*asinh(If0./(2*i0n));
    posEta3 = 2*R*T/F*asinh(If3./(2*i0p));

    % Compute cell voltage (ROMout.Vcell)
    Uocpn0 = cellData.function.neg.Uocp(Z(ind.Thetass0),T);
    Uocpp3 = cellData.function.pos.Uocp(Z(ind.Thetass3),T);
    Rfn    = cellData.function.neg.Rf(SOCnAvg,T);
    Rfp    = cellData.function.pos.Rf(SOCpAvg,T);

    Vcell = posEta3 - negEta0 + PhieTilde3 + Uocpp3 - Uocpn0 ...
            + (Rfp*Z(ind.Ifdl3) - Rfn*Z(ind.Ifdl0));

    % Compute solid potential (phis)
    Z(ind.posPhis) = Z(ind.posPhis) + Vcell;
    
    % Finally, compute SOC 
    Zsoc = ekfData.SOC0 - x0*(ekfData.Ts/(3600*ekfData.Q));     
  end

  %% ----------------------------------------------------------------------
  % This function computes derivative of voltage with respect to states
  function [Chat,Chat0] = getChatV(Xind,Z,Tk)
    % Full-size "C" matrices from the four nearest-neighbor models
    C1 = Xind.gamma(1)*ekfData.M(Xind.theT(1),Xind.theZ(1)).C;
    C2 = Xind.gamma(2)*ekfData.M(Xind.theT(2),Xind.theZ(2)).C;
    C3 = Xind.gamma(3)*ekfData.M(Xind.theT(3),Xind.theZ(3)).C;
    C4 = Xind.gamma(4)*ekfData.M(Xind.theT(4),Xind.theZ(4)).C;

    % Voltage has several parts: Each is addressed below by number
    % 1) Rfp * Ifdlp - Rfn * Ifdln
    % 2) etap - etan, approximated as Rctp * Ifp - Rctn * Ifn
    % 3) phietildep
    % 4) Uocpp(thetassp) - Uocpn(thetassn) [only term with integrator...]
    
    % 1) Terms for film resistance and Ifdl
    switch ekfData.method
      case 'OB'
        x0 = ekfData.x0;
        xSOC = ekfData.SOC0 - x0*(ekfData.Ts/(3600*ekfData.Q)); 
        SOCnAvg = cellData.function.neg.soc(xSOC,Tk);
        SOCpAvg = cellData.function.pos.soc(xSOC,Tk);
        Rfn = cellData.function.neg.Rf(SOCnAvg,Tk);
        Rfp = cellData.function.pos.Rf(SOCpAvg,Tk);
        
        Chat{1} = Rfp*C1(ind.Ifdl3,:) - Rfn*C1(ind.Ifdl0,:);
        Chat{2} = Rfp*C2(ind.Ifdl3,:) - Rfn*C2(ind.Ifdl0,:);
        Chat{3} = Rfp*C3(ind.Ifdl3,:) - Rfn*C3(ind.Ifdl0,:);
        Chat{4} = Rfp*C4(ind.Ifdl3,:) - Rfn*C4(ind.Ifdl0,:);
      case 'MB'
        x0 = ekfData.xhat(end);
        xSOC = ekfData.SOC0 - x0*(ekfData.Ts/(3600*ekfData.Q)); 
        SOCnAvg = cellData.function.neg.soc(xSOC,Tk);
        SOCpAvg = cellData.function.pos.soc(xSOC,Tk);
        Rfn = cellData.function.neg.Rf(SOCnAvg,Tk);
        Rfp = cellData.function.pos.Rf(SOCpAvg,Tk);
        
        Chat = Rfp*sum([C1(ind.Ifdl3,:); C2(ind.Ifdl3,:); ...
                        C3(ind.Ifdl3,:); C4(ind.Ifdl3,:)]) ...
                - Rfn*sum([C1(ind.Ifdl0,:); C2(ind.Ifdl0,:); ...
                           C3(ind.Ifdl0,:); C4(ind.Ifdl0,:)]);
    end
    
    % 2) Terms for charge-transfer resistance and eta (simplified method)
    k0n = cellData.function.neg.k0(SOCnAvg,Tk);
    k0p = cellData.function.pos.k0(SOCpAvg,Tk);
    i0n = k0n*sqrt(Z(ind.Thetae(1))*(1-Z(ind.Thetass0))*Z(ind.Thetass0));
    i0p = k0p*sqrt(Z(ind.Thetae(end))*(1-Z(ind.Thetass3))*Z(ind.Thetass3));
    Rctn = R*Tk/(F*i0n);
    Rctp = R*Tk/(F*i0p);    
    switch ekfData.method
      case 'OB'
        Chat{1} = Chat{1} + Rctp*C1(ind.If3,:) - Rctn*C1(ind.If0,:);
        Chat{2} = Chat{2} + Rctp*C2(ind.If3,:) - Rctn*C2(ind.If0,:);
        Chat{3} = Chat{3} + Rctp*C3(ind.If3,:) - Rctn*C3(ind.If0,:);
        Chat{4} = Chat{4} + Rctp*C4(ind.If3,:) - Rctn*C4(ind.If0,:);
      case 'MB'
        Chat = Chat + Rctp*sum([C1(ind.If3,:); C2(ind.If3,:); 
                                C3(ind.If3,:); C4(ind.If3,:)]) ...
                - Rctn*sum([C1(ind.If0,:); C2(ind.If0,:); 
                            C3(ind.If0,:); C4(ind.If0,:)]);
    end
    
    % 3) Terms for debiased electrolyte potential drop
    switch ekfData.method
      case 'OB'
        Chat{1} = Chat{1} + C1(ind.Phie(end),:);
        Chat{2} = Chat{2} + C2(ind.Phie(end),:);
        Chat{3} = Chat{3} + C3(ind.Phie(end),:);
        Chat{4} = Chat{4} + C4(ind.Phie(end),:);
      case 'MB'
        Chat = Chat + sum([C1(ind.Phie(end),:); C2(ind.Phie(end),:); 
                           C3(ind.Phie(end),:); C4(ind.Phie(end),:)]);
    end
    
    % 4) Terms for Uocp
    dUocpn0 = cellData.function.neg.dUocp(Z(ind.Thetass0),Tk);
    dUocpp3 = cellData.function.pos.dUocp(Z(ind.Thetass3),Tk);
    res0n   = -dUocpn0*ekfData.Ts*(cellData.function.neg.soc(1,Tk) - ...
               cellData.function.neg.soc(0,Tk))/(3600*ekfData.Q); 
    res0p   = -dUocpp3*ekfData.Ts*(cellData.function.pos.soc(1,Tk) - ...
               cellData.function.pos.soc(0,Tk))/(3600*ekfData.Q);    
    Chat0 = res0p - res0n;
    switch ekfData.method
      case 'OB'
        Chat{1} = Chat{1} + (dUocpp3*C1(ind.Thetass3,:) ...
                      - dUocpn0*C1(ind.Thetass0,:));
        Chat{2} = Chat{2} + (dUocpp3*C2(ind.Thetass3,:) ...
                      - dUocpn0*C2(ind.Thetass0,:));
        Chat{3} = Chat{3} + (dUocpp3*C3(ind.Thetass3,:) ...
                      - dUocpn0*C3(ind.Thetass0,:));
        Chat{4} = Chat{4} + (dUocpp3*C4(ind.Thetass3,:) ...
                      - dUocpn0*C4(ind.Thetass0,:));
      case 'MB'
        Chat = Chat + dUocpp3*sum([C1(ind.Thetass3,:); C2(ind.Thetass3,:); 
                              C3(ind.Thetass3,:); C4(ind.Thetass3,:)]) ...
                - dUocpn0*sum([C1(ind.Thetass0,:); C2(ind.Thetass0,:); 
                              C3(ind.Thetass0,:); C4(ind.Thetass0,:)]);
        Chat = [Chat Chat0]; % as final step, augment integrator state
    end
  end

  %% ----------------------------------------------------------------------
  % This function computes derivative of voltage with respect to states
  function [Chat,ChatV,Chat0,ChatV0] = getChatZ(Xind,Z,Tk)
    % Full-size "C" matrices from the four nearest-neighbor models
    C1 = Xind.gamma(1)*ekfData.M(Xind.theT(1),Xind.theZ(1)).C;
    C2 = Xind.gamma(2)*ekfData.M(Xind.theT(2),Xind.theZ(2)).C;
    C3 = Xind.gamma(3)*ekfData.M(Xind.theT(3),Xind.theZ(3)).C;
    C4 = Xind.gamma(4)*ekfData.M(Xind.theT(4),Xind.theZ(4)).C;

    % Different variables need different modifications to model "C" matrix
    
    % Variables Ifdl, If, Thetae and Phisn do not require changes to "C"
    % Simply copy "C" matrices over
    switch ekfData.method
      case 'OB'
        Chat{1} = C1; Chat{2} = C2; Chat{3} = C3; Chat{4} = C4;
      case 'MB'
        Chat = C1 + C2 + C3 + C4;
    end
    Chat0 = zeros(size(C1,1),1); % initialize to no residue terms

    % Phis in the positive electrode has a "vcell" component
    [ChatV,ChatV0] = getChatV(Xind,Z,Tk);
    switch ekfData.method
      case 'OB'
        for kk = 1:4
          C = Chat{kk};
          for thePhis = 1:length(ind.posPhis)
            C(ind.posPhis(thePhis),:) = C(ind.posPhis(thePhis),:) + ChatV{kk};
          end
          Chat{kk} = C;
        end
        Chat0(ind.posPhis) = ChatV0;
      case 'MB'
        for thePhis = 1:length(ind.posPhis)
          Chat(ind.posPhis(thePhis),:) = Chat(ind.posPhis(thePhis),:) + ChatV(1:end-1);
        end
        Chat0(ind.posPhis) = ChatV(end);
    end
    
    % Thetass includes a term that affects "C" for x0.
    res0n   = -ekfData.Ts*(cellData.function.neg.soc(1,Tk) - ...
               cellData.function.neg.soc(0,Tk))/(3600*ekfData.Q); 
    res0p   = -ekfData.Ts*(cellData.function.pos.soc(1,Tk) - ...
               cellData.function.pos.soc(0,Tk))/(3600*ekfData.Q);    
    Chat0(ind.negThetass) = res0n;
    Chat0(ind.posThetass) = res0p;
    
    % Phise includes a term that affects "C" for x0.
    switch ekfData.method
      case 'OB'
        x0 = ekfData.x0;
      case 'MB'
        x0 = ekfData.xhat(end);
    end
    xSOC = ekfData.SOC0 - x0*(ekfData.Ts/(3600*ekfData.Q)); 
    SOCnAvg = cellData.function.neg.soc(xSOC,Tk);
    SOCpAvg = cellData.function.pos.soc(xSOC,Tk);
    dUocpn  = cellData.function.neg.dUocp(SOCnAvg,Tk);
    dUocpp  = cellData.function.pos.dUocp(SOCpAvg,Tk);
    
    Chat0(ind.negPhise) = dUocpn*res0n;
    Chat0(ind.posPhise) = dUocpp*res0p;
    
    % Phie includes Phise(0,t) term and also integrator term
    Chat0(ind.Phie) = -dUocpn*res0n;
    switch ekfData.method
      case 'OB'
        for kk = 1:4
          C = Chat{kk};
          for thePhie = 1:length(ind.Phie)
            C(ind.Phie(thePhie),:) = C(ind.Phie(thePhie),:) - C(ind.Phise0,:);
          end
          Chat{kk} = C;
        end
      case 'MB'
        for thePhie = 1:length(ind.Phie)
          Chat(ind.Phie(thePhie),:) = Chat(ind.Phie(thePhie),:) - Chat(ind.Phise0,:);
        end        
        Chat = [Chat, Chat0]; % as final step, augment integrator state
    end
  end

  %% ----------------------------------------------------------------------
  % This function sets up indices (ind) into the model linear output
  % vector "y" of variables needed to compute cell voltage and nonlinear
  % corrections. It also sets up locations (loc) in xtilde coordinates
  % of electrolyte variables.
  % -----------------------------------------------------------------------
  function setupIndsLocs
    % -- Find indices of outputs in model structure, to be used later 
    tfName = ekfData.tfData.names; % TF names
    tfLocs = ekfData.tfData.xLoc;  % TF regions (normalized x locations)

    % Negative electrode
    ind.negIfdl    = find(strcmp(tfName,'negIfdl') == 1);
    ind.negIf      = find(strcmp(tfName,'negIf') == 1);
    ind.negPhis    = find(strcmp(tfName,'negPhis') == 1);
    ind.negPhise   = find(strcmp(tfName,'negPhise') == 1);
    ind.negThetass = find(strcmp(tfName,'negThetass') == 1);
    loc.negIfdl    = tfLocs(ind.negIfdl);
    loc.negIf      = tfLocs(ind.negIf);
    loc.negPhis    = tfLocs(ind.negPhis);
    loc.negPhise   = tfLocs(ind.negPhise);
    loc.negThetass = tfLocs(ind.negThetass);

    % Positive electrode
    ind.posIfdl    = find(strcmp(tfName,'posIfdl') == 1);
    ind.posIf      = find(strcmp(tfName,'posIf') == 1);
    ind.posPhis    = find(strcmp(tfName,'posPhis') == 1);
    ind.posPhise   = find(strcmp(tfName,'posPhise') == 1);
    ind.posThetass = find(strcmp(tfName,'posThetass') == 1);
    loc.posIfdl    = tfLocs(ind.posIfdl);
    loc.posIf      = tfLocs(ind.posIf);
    loc.posPhis    = tfLocs(ind.posPhis);
    loc.posPhise   = tfLocs(ind.posPhise);
    loc.posThetass = tfLocs(ind.posThetass);

    % Electrolyte potential across cell width
    ind.negPhie = find(strcmp(tfName,'negPhie') == 1);
    ind.sepPhie = find(strcmp(tfName,'sepPhie') == 1); 
    ind.posPhie = find(strcmp(tfName,'posPhie') == 1);
    loc.negPhie = tfLocs(ind.negPhie);
    loc.sepPhie = tfLocs(ind.sepPhie); 
    loc.posPhie = tfLocs(ind.posPhie);

    % Electrolyte normalized concentration across cell width
    ind.negThetae  = find(strcmp(tfName,'negThetae')== 1);
    ind.sepThetae  = find(strcmp(tfName,'sepThetae')== 1); 
    ind.posThetae  = find(strcmp(tfName,'posThetae')== 1);
    loc.negThetaes = tfLocs(ind.negThetae);
    loc.sepThetaes = tfLocs(ind.sepThetae); 
    loc.posThetaes = tfLocs(ind.posThetae);

    % Combine inds and locs for variables across entire cell width
    ind.Ifdl    = [ind.negIfdl; ind.posIfdl];
    loc.Ifdl    = [loc.negIfdl; loc.posIfdl];
    ind.If      = [ind.negIf; ind.posIf];
    loc.If      = [loc.negIf; loc.posIf];
    ind.Phis    = [ind.negPhis; ind.posPhis];
    loc.Phis    = [loc.negPhis; loc.posPhis];
    ind.Phise   = [ind.negPhise; ind.posPhise];
    loc.Phise   = [loc.negPhise; loc.posPhise];
    ind.Thetass = [ind.negThetass;ind.posThetass];
    loc.Thetass = [loc.negThetass;loc.posThetass];
    
    ind.Phie    = [ind.negPhie;ind.sepPhie;ind.posPhie];
    loc.Phie    = [loc.negPhie;loc.sepPhie;loc.posPhie];
    ind.Thetae  = [ind.negThetae;ind.sepThetae;ind.posThetae];
    loc.Thetae  = [loc.negThetaes;loc.sepThetaes;loc.posThetaes];

    %-- Verify that variables required for nonlinear corrections exist 
    % Need to check:
    %  1. Ifdl    at both current-collectors
    %  2. If   at both current-collectors
    %  3. ROMout.Thetae  at both current-collectors
    %  4. Thetass at both current-collectors
    %  5. Phise   at negative-electrode current-collector 
    %  6. ROMout.Phie    at positive-electrode current-collector

    % Find location indexes
    ind.Ifdl0    = ind.Ifdl(loc.Ifdl == 0);
    ind.Ifdl3    = ind.Ifdl(loc.Ifdl == 3);
    ind.If0      = ind.If(loc.If == 0);
    ind.If3      = ind.If(loc.If == 3);
    ind.Thetass0 = ind.Thetass(loc.Thetass == 0);
    ind.Thetass3 = ind.Thetass(loc.Thetass == 3);
    ind.Phise0   = ind.Phise(loc.Phise == 0);
    ind.Thetae0  = ind.Thetae(loc.Thetae == 0);
    ind.Thetae3  = ind.Thetae(abs(loc.Thetae - 3) < 0.001);

    % Check #1
    if isempty(ind.Ifdl0)
      error('Simulation requires ifdl at negative-collector!'); 
    end
    if isempty(ind.Ifdl3)
      error('Simulation requires ifdl at positive-collector!'); 
    end

    % Check #2
    if isempty(ind.If0)
      error('Simulation requires if at negative-collector!'); 
    end
    if isempty(ind.If3)
      error('Simulation requires if at positive-collector!'); 
    end

    % Check #3
    if loc.Thetae(1) > 0
      error('Simulation requires thetae at negative-collector!'); end
    if or(loc.Thetae(end)>3+eps,loc.Thetae(end)<3-eps)
      error('Simulation requires thetae at positive-collector!'); end

    % Check #4
    if isempty(ind.Thetass0)
      error('Simulation requires thetass at negative-collector!'); 
    end
    if isempty(ind.Thetass3)
      error('Simulation requires thetass at positive-collector!'); 
    end

    % Check #5
    if isempty(ind.Phise0)
      error('Simulation requires phise at negative-collector!'); 
    end

    % Check #6
    if loc.Phie(1) == 0
      shortWarn('First phie x-location should not be zero. Ignoring');
      ind.Phie = ind.Phie(2:end); loc.Phie = loc.Phie(2:end);
    end
    if or(loc.Phie(end)>3+eps,loc.Phie(end)<3-eps)
      error('Simulation requires phie at positive-collector!');
    end
  end
 
  %% ----------------------------------------------------------------------
  % This function displays a short warning message to the user if warnState
  % is []. (Does not display line numbers of warning, etc.)
  % -----------------------------------------------------------------------
  function shortWarn(msg)
    persistent warnState
    if strcmpi(msg,'on')
      warnState = []; 
    elseif strcmpi(msg,'off')
      warnState = 1;
    elseif isempty(warnState)
      cprintf([1,1/2,0],[' - Warning: ' msg '\n']);
      warnCount = warnCount + 1;
    end
  end

end