% function ROMout = outBlendPack(ROMs,simData,architecture)
% 
% Inputs:
%   ROMs         = matrix of reduced-order models created with an XRA
%   simData      = simulation profile loaded using "loadInput"
%   architecture = either 'SCM' or 'PCM'; the pack configuraion
%
% Outputs:
%   ROMout       = output structure contains simulated info
%
% This function simulates an Ns by Np based battery pack using the
% physics-based ROM and output blending 
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

function ROMout = outBlendPack(ROMs,simData,architecture)
  shortWarn('off');
  
  % Re-define current and temperature data to enforce fixed Ts 
  [tk,ik,~,Tk] = resampleInput;
  
  [Ns,Np] = size(ROMs.ROM); % size of battery pack
  
  % Set up indices of simulation variables in linear outputs "y" 
  setupIndsLocs; % stores in ROMs.ind and ROMs.loc
  
  % Create storage for simulation variables in output data structure
  setupOutputs;  % ROMout(1:Ns,1:Np).(fields)
  
  for ks = 1:Ns
    for kp = 1:Np
      % Set up initial electrode local SOC (convert from percent)
      ROMs.ROM(ks,kp).SOC0n = ROMs.ROM(ks,kp).cellData.function.neg.soc(...
        simData.SOC0(ks,kp)/100,Tk(1));
      ROMs.ROM(ks,kp).SOC0p = ROMs.ROM(ks,kp).cellData.function.pos.soc(...
        simData.SOC0(ks,kp)/100,Tk(1));

      % Set up blending method and associated matrices/vectors
      [bigA,bigX,Tspts,Zspts,ZZ] = setupBlend(ks,kp);
      ROMs.ROM(ks,kp).bigA = bigA;
      ROMs.ROM(ks,kp).Tspts = Tspts;
      ROMs.ROM(ks,kp).Zspts = Zspts;
      ROMs.ROM(ks,kp).ZZ = ZZ;

      % Set up the cell state
      ROMs.ROM(ks,kp).cellState.bigX = bigX;
      ROMs.ROM(ks,kp).cellState.SOCnAvg = ROMs.ROM(ks,kp).SOC0n;
      ROMs.ROM(ks,kp).cellState.SOCpAvg = ROMs.ROM(ks,kp).SOC0p;
    end
  end
  
  options = optimset('Display','off');
  
  % Main simulation loop
  fprintf('Simulating ROM via outBlend...\n'); 
  for k = 0:length(ik)-1
    switch upper(architecture)
      case 'PCM'
        for ks = 1:Ns
          ij0 = ik(k+1)/Np*ones(Np-1,1);
          ijk = fminsearch(@optPCM,ij0,options);
          updatePCM(ijk);
        end        
      case 'SCM'
        ij0 = ik(k+1)/Np*ones(Np-1,1);
        ijk = fminsearch(@optSCM,ij0,options);
        updateSCM(ijk);
      otherwise
        error('Unknown architecture %s to outBlendPack.m',architecture)
    end
    % Update information to user every 100 iterations
    if (rem(k,100) == 0)
      fprintf('%d %2.8f\n',k,ROMout(1,1).cellSOC(k+1)); 
    end    
  end
  
  % Return values as a structure
  storeData;

  % Finished! Now, return to the user

  %% ======================================================================
  % The functions below this point implement the details of the higher-
  % level functionality indicated above
  % =======================================================================

  %% ----------------------------------------------------------------------
  % This function resamples the simulation time/temperature/current vectors
  % to ensure an even sampling interval.
  % -----------------------------------------------------------------------
  function [tk,ik,pk,Tk] = resampleInput
    Ts = ROMs.ROM(1).xraData.Tsamp; 
    tk = simData.time(1):Ts:simData.time(end); 
      ik = interp1(simData.time,simData.Iapp,tk); 
      pk = 0*ik;
    Tk = interp1(simData.time,simData.T,tk)+273.15; 
    if size(ik,1) ~= size(Tk,1)
      error('Current and temperature profiles must be same length');
    end
  end

  %% ----------------------------------------------------------------------
  % This function sets up indices (ROMs.ind) into the model linear output
  % vector "y" of variables needed to compute cell voltage and nonlinear
  % corrections. It also sets up locations (ROMs.loc) in xtilde coordinates
  % of electrolyte variables.
  % -----------------------------------------------------------------------
  function setupIndsLocs
    % -- Find indices of outputs in model structure, to be used later 
    tfName = ROMs.ROM(1,1).tfData.names; % TF names
    tfLocs = ROMs.ROM(1,1).tfData.xLoc;  % TF regions (normalized x locations)
    ROMs.tfLocs = tfLocs;

    % Negative electrode
    ROMs.ind.negIfdl    = find(strcmp(tfName,'negIfdl') == 1);
    ROMs.ind.negIf      = find(strcmp(tfName,'negIf') == 1);
    ROMs.ind.negIdl     = find(strcmp(tfName,'negIdl') == 1);
    ROMs.ind.negPhis    = find(strcmp(tfName,'negPhis') == 1);
    ROMs.ind.negPhise   = find(strcmp(tfName,'negPhise') == 1);
    ROMs.ind.negThetass = find(strcmp(tfName,'negThetass') == 1);

    % Positive electrode
    ROMs.ind.posIfdl    = find(strcmp(tfName,'posIfdl') == 1);
    ROMs.ind.posIf      = find(strcmp(tfName,'posIf') == 1);
    ROMs.ind.posIdl     = find(strcmp(tfName,'posIdl') == 1);
    ROMs.ind.posPhis    = find(strcmp(tfName,'posPhis') == 1);
    ROMs.ind.posPhise   = find(strcmp(tfName,'posPhise') == 1);
    ROMs.ind.posThetass = find(strcmp(tfName,'posThetass') == 1);

    % Electrolyte potential across cell width
    ROMs.ind.negPhie = find(strcmp(tfName,'negPhie') == 1);
    ROMs.ind.sepPhie = find(strcmp(tfName,'sepPhie') == 1); 
    ROMs.ind.posPhie = find(strcmp(tfName,'posPhie') == 1);
    ROMs.loc.negPhie = tfLocs(ROMs.ind.negPhie);
    ROMs.loc.sepPhie = tfLocs(ROMs.ind.sepPhie); 
    ROMs.loc.posPhie = tfLocs(ROMs.ind.posPhie);

    ROMs.ind.Phie = [ROMs.ind.negPhie;ROMs.ind.sepPhie;ROMs.ind.posPhie];
    ROMs.loc.Phie = [ROMs.loc.negPhie;ROMs.loc.sepPhie;ROMs.loc.posPhie];

    % Electrolyte normalized concentration across cell width
    ROMs.ind.negThetae  = find(strcmp(tfName,'negThetae')== 1);
    ROMs.ind.sepThetae  = find(strcmp(tfName,'sepThetae')== 1); 
    ROMs.ind.posThetae  = find(strcmp(tfName,'posThetae')== 1);
    ROMs.loc.negThetaes = tfLocs(ROMs.ind.negThetae);
    ROMs.loc.sepThetaes = tfLocs(ROMs.ind.sepThetae); 
    ROMs.loc.posThetaes = tfLocs(ROMs.ind.posThetae);

    ROMs.ind.Thetae = [ROMs.ind.negThetae;ROMs.ind.sepThetae;...
                       ROMs.ind.posThetae];
    ROMs.loc.Thetae = [ROMs.loc.negThetaes;ROMs.loc.sepThetaes;...
                       ROMs.loc.posThetaes];

    %-- Verify that variables required for nonlinear corrections exist 
    % Need to check:
    %  1. Ifdl    at both current-collectors
    %  2. If   at both current-collectors
    %  3. ROMout.Thetae  at both current-collectors
    %  4. Thetass at both current-collectors
    %  5. Phise   at negative-electrode current-collector 
    %  6. ROMout.Phie    at positive-electrode current-collector

    % Find location indexes
    ROMs.ind.negIfdl0    = find(strcmp(tfName,'negIfdl') == 1 ...
                                & tfLocs == 0);
    ROMs.ind.posIfdl3    = find(strcmp(tfName,'posIfdl') == 1 ...
                                & tfLocs == 3);
    ROMs.ind.negIf0   = find(strcmp(tfName,'negIf') == 1 ...
                                & tfLocs == 0);
    ROMs.ind.posIf3   = find(strcmp(tfName,'posIf') == 1 ...
                                & tfLocs == 3);
    ROMs.ind.negThetass0 = find(strcmp(tfName,'negThetass') == 1 ...
                                & tfLocs == 0);
    ROMs.ind.posThetass3 = find(strcmp(tfName,'posThetass') == 1 ...
                                & tfLocs == 3);
    ROMs.ind.negPhise0   = find(strcmp(tfName,'negPhise') == 1 ...
                                & tfLocs == 0);

    % Check #1
    if isempty(ROMs.ind.negIfdl0)
      error('Simulation requires ifdl at negative-collector!'); 
    end
    if isempty(ROMs.ind.posIfdl3)
      error('Simulation requires ifdl at positive-collector!'); 
    end

    % Check #2
    if isempty(ROMs.ind.negIf0)
      error('Simulation requires if at negative-collector!'); 
    end
    if isempty(ROMs.ind.posIf3)
      error('Simulation requires if at positive-collector!'); 
    end

    % Check #3
    if ROMs.loc.Thetae(1) > 0
      error('Simulation requires thetae at negative-collector!'); end
    if or(ROMs.loc.Thetae(end)>3+eps,ROMs.loc.Thetae(end)<3-eps)
      error('Simulation requires thetae at positive-collector!'); end

    % Check #4
    if isempty(ROMs.ind.negThetass0)
      error('Simulation requires thetass at negative-collector!'); 
    end
    if isempty(ROMs.ind.posThetass3)
      error('Simulation requires thetass at positive-collector!'); 
    end

    % Check #5
    if isempty(ROMs.ind.negPhise0)
      error('Simulation requires phise at negative-collector!'); 
    end

    % Check #6
    if ROMs.loc.Phie(1) == 0
      shortWarn('First phie x-location should not be zero. Ignoring');
      ROMs.ind.Phie = ROMs.ind.Phie(2:end);
    end
    if or(ROMs.loc.Phie(end)>3+eps,ROMs.loc.Phie(end)<3-eps)
      error('Simulation requires phie at positive-collector!');
    end
  end

  %% ----------------------------------------------------------------------
  % This function initializes storage for simulation variables in the
  % output data structure
  % -----------------------------------------------------------------------
  function setupOutputs
    duration = length(ik);
    
    for kks = 1:Ns
      for kkp = 1:Np        
        % At negative electrode and its current collector
        ROMout(kks,kkp).negIfdl    = zeros(duration,length(ROMs.ind.negIfdl));
        ROMout(kks,kkp).negIf      = zeros(duration,length(ROMs.ind.negIf));
        ROMout(kks,kkp).negIdl     = zeros(duration,length(ROMs.ind.negIdl));
        ROMout(kks,kkp).negPhis    = zeros(duration,length(ROMs.ind.negPhis));
        ROMout(kks,kkp).negPhise   = zeros(duration,length(ROMs.ind.negPhise));
        ROMout(kks,kkp).negThetass = zeros(duration,length(ROMs.ind.negThetass));

        ROMout(kks,kkp).negIfdl0    = zeros(duration,length(ROMs.ind.negIfdl0));
        ROMout(kks,kkp).negIf0      = zeros(duration,length(ROMs.ind.negIf0));
        ROMout(kks,kkp).negEta0     = zeros(duration,length(ROMs.ind.negIf0));
        ROMout(kks,kkp).negPhise0   = zeros(duration,length(ROMs.ind.negPhise0));
        ROMout(kks,kkp).negThetass0 = zeros(duration,length(ROMs.ind.negThetass0));

        % At positive electrode and its current collector
        ROMout(kks,kkp).posIfdl    = zeros(duration,length(ROMs.ind.posIfdl));
        ROMout(kks,kkp).posIf      = zeros(duration,length(ROMs.ind.posIf));
        ROMout(kks,kkp).posIdl     = zeros(duration,length(ROMs.ind.posIdl));
        ROMout(kks,kkp).posPhis    = zeros(duration,length(ROMs.ind.posPhis));
        ROMout(kks,kkp).posPhise   = zeros(duration,length(ROMs.ind.posPhise));
        ROMout(kks,kkp).posThetass = zeros(duration,length(ROMs.ind.posThetass));

        ROMout(kks,kkp).posIfdl3    = zeros(duration,length(ROMs.ind.posIfdl3));
        ROMout(kks,kkp).posIf3      = zeros(duration,length(ROMs.ind.posIf3));
        ROMout(kks,kkp).posEta3     = zeros(duration,length(ROMs.ind.posIf3));
        ROMout(kks,kkp).posThetass3 = zeros(duration,length(ROMs.ind.posThetass3));

        % Electrolyte variables
        ROMout(kks,kkp).Phie   = zeros(duration,length(ROMs.ind.Phie)+1); % add x=0 loc
        ROMout(kks,kkp).Thetae = zeros(duration,length(ROMs.ind.Thetae));

        % Other cell and electrode quantities
        ROMout(kks,kkp).Vcell   = zeros(duration,1);
        ROMout(kks,kkp).cellSOC = zeros(duration,1); 
        ROMout(kks,kkp).negSOC  = zeros(duration,1); 
        ROMout(kks,kkp).posSOC  = zeros(duration,1); 
      end
    end
  end      

  %% ----------------------------------------------------------------------
  % This function sets up the internal structure for performing the output
  % blending for the cell at (kks,kkp)
  % -----------------------------------------------------------------------
  function [bigA,bigX,Tspts,Zspts,ZZ] = setupBlend(kks,kkp)
    % Pluck out the ROM to work with
    ROM = ROMs.ROM(kks,kkp);
    
    % Extra all temperature and SOC setpoints (sort in ascending order)
    Tspts = sort(ROM.xraData.T)+273.15;
    Zspts = sort(ROM.xraData.SOC/100);

    % Create a "bigA" matrix.
    % Each column is the diagonal values of a specific A matrix.
    % The setpoint convention from column 1 to the last column is:
    % Z(1)T(1), Z(2)T(1), ..., Z(1)T(2), Z(2)T(2),..., Z(3)T(1),...
    % ROMmdls generated by xRA is a #Tspts by # Zspts structure
    [TT,ZZ] = size(ROM.ROMmdls);
    bigA = zeros(length(ROM.ROMmdls(1,1).A),TT*ZZ); 
    for tt = 1:TT
      for zz = 1:ZZ
        ind  = (tt-1)*ZZ + zz;
        ROMi = ROM.ROMmdls(tt,zz);
        bigA(:,ind) = real(diag(ROMi.A));

        % Don't forget to change res0 because we need only [phise]*
        ROMs.ROM(kks,kkp).ROMmdls(tt,zz).C(ROMs.ind.negPhise,end) = 0;
        ROMs.ROM(kks,kkp).ROMmdls(tt,zz).C(ROMs.ind.posPhise,end) = 0;
      end
    end

    % Create the "bigX" matrix.
    % Each column stores the system states at a specific setpoint.
    bigX = zeros(size(bigA));
  end

  %% ----------------------------------------------------------------------
  % This function produces a cost value that must be minimized during the
  % simulation at every timestep to ensure that cell currents in a PCM add
  % up to pack current and all cell voltages in a PCM are identical
  % -----------------------------------------------------------------------
  function cost = optPCM(ij0)
    ijvec = [ij0(:); ik(k+1) - sum(ij0(:))]; % enforce total current
    Vcell = 0*ijvec;
    for kkp = 1:Np
      Vcell(kkp) = simStep(ijvec(kkp),Tk(k+1),ks,kkp);
    end
    cost = rms(diff(Vcell)); % small cost if voltages are similar
  end

  %% ----------------------------------------------------------------------
  % This function updates all cells in a PCM for a given vector of
  % currents, enforcing that the sum of currents equals the pack current
  % -----------------------------------------------------------------------
  function updatePCM(ij0)
    ijvec = [ij0(:); ik(k+1) - sum(ij0(:))]; % enforce total current
    for kkp = 1:Np
      [~,ROMs.ROM(ks,kkp).cellState] = simStep(ijvec(kkp),Tk(k+1),ks,kkp);
    end
  end

  %% ----------------------------------------------------------------------
  % This function produces a cost value that must be minimized during the
  % simulation at every timestep to ensure that SCM currents add up to pack
  % current and all SCM total voltages are identical 
  % -----------------------------------------------------------------------
  function cost = optSCM(ij0)
    ijvec = [ij0(:); ik(k+1) - sum(ij0(:))]; % enforce total current
    Vscm = 0*ijvec;
    for kks = 1:Ns
      for kkp = 1:Np
        Vscm(kkp) = Vscm(kkp) + simStep(ijvec(kkp),Tk(k+1),kks,kkp);
      end
    end
    cost = rms(diff(Vscm)); % small cost if voltages are similar
  end

  %% ----------------------------------------------------------------------
  % This function updates all cells in a SCM-based pack for a given vector
  % of SCM currents, enforcing that the sum of currents equals the pack
  % current 
  % -----------------------------------------------------------------------
  function updateSCM(ij0)
    ijvec = [ij0(:); ik(k+1) - sum(ij0(:))]; % enforce total current
    for kks = 1:Ns
      for kkp = 1:Np
        [~,ROMs.ROM(kks,kkp).cellState] = simStep(ijvec(kkp),Tk(k+1),kks,kkp);
      end
    end
  end

  %% ----------------------------------------------------------------------
  % This function simulates one time step using output blending. Note that
  % while data are stored in the (k+1)st index, they can be overwritten if
  % the calling routine decides that an over/under current/power condition
  % has occured and re-calls this function
  % -----------------------------------------------------------------------
  function [Vcell,newCellState] = simStep(Iapp,T,kks,kkp)
    theROM    = ROMs.ROM(kks,kkp);
    cellState = theROM.cellState;
    cellData  = theROM.cellData;
    ROMmdls   = theROM.ROMmdls;
    SOC0n     = theROM.SOC0n;
    SOC0p     = theROM.SOC0p;
    
    F         = cellData.const.F;
    R         = cellData.const.R;    
    Q         = cellData.function.const.Q(); % cell capacity in Ah
    Rc        = cellData.function.const.Rc();
    theta0n   = cellData.function.neg.theta0();
    theta0p   = cellData.function.pos.theta0();
    theta100n = cellData.function.neg.theta100();
    theta100p = cellData.function.pos.theta100();
    
    wDLn      = cellData.function.neg.wDL(SOC0n,T); 
    wDLp      = cellData.function.pos.wDL(SOC0p,T); 
    Cdln      = cellData.function.neg.Cdl(SOC0n,T); 
    Cdlp      = cellData.function.pos.Cdl(SOC0p,T);
    nDLn      = cellData.function.neg.nDL();
    nDLp      = cellData.function.pos.nDL();
    Cdleffn   = (Cdln^(2-nDLn))*(wDLn^(nDLn-1));
    Cdleffp   = (Cdlp^(2-nDLp))*(wDLp^(nDLp-1));

    % Load present model state from "cellState"
    bigX      = cellState.bigX;
    SOCnAvg   = cellState.SOCnAvg;
    SOCpAvg   = cellState.SOCpAvg;
    
    % ------------------------------------------
    % Step 1: Update cell SOC for next time step
    % ------------------------------------------
    % Store electrode average SOC data and compute cell SOC
    ROMout(kks,kkp).negSOC(k+1)  = SOCnAvg;
    ROMout(kks,kkp).posSOC(k+1)  = SOCpAvg;
    ROMout(kks,kkp).cellSOC(k+1) = (SOCnAvg - theta0n)/(theta100n - theta0n);

    % Compute integrator input gains
    dUocpnAvg = cellData.function.neg.dUocp(SOCnAvg,T); 
    dUocppAvg = cellData.function.pos.dUocp(SOCpAvg,T);  
    dQn       = abs(theta100n-theta0n);
    dQp       = abs(theta100p-theta0p);   
    res0n     = -dQn/(3600*Q-Cdleffn*dQn*dUocpnAvg);
    res0p     =  dQp/(3600*Q-Cdleffp*dQp*dUocppAvg);
    
    % Compute average negative-electrode SOC
    SOCnAvg = SOCnAvg + res0n*Iapp*(tk(2)-tk(1));
    if SOCnAvg < 0, shortWarn('Average SOCn < 0'); SOCnAvg = 0; end
    if SOCnAvg > 1, shortWarn('Average SOCn > 1'); SOCnAvg = 1; end

    % Compute average positive-electrode SOC
    SOCpAvg = SOCpAvg + res0p*Iapp*(tk(2)-tk(1));
    if SOCpAvg < 0, shortWarn('Average SOCp < 0'); SOCpAvg = 0; end
    if SOCpAvg > 1, shortWarn('Average SOCp > 1'); SOCpAvg = 1; end

    % ---------------------------------------------------------------------
    % Step 2: Update linear output of closest four pre-computed ROM as
    %         y[k] = C*x[k] + D*iapp[k]
    % ---------------------------------------------------------------------
    % Find the two closest Zspts setpoints: "Zupper" and "Zlower"
    dZ = abs(ROMout(kks,kkp).cellSOC(k+1)-Zspts); [~,iZ] = sort(dZ);
    Zupper = Zspts; iZupper = iZ;
    Zlower = Zspts; iZlower = iZ;
    if length(iZ)>1
      Zupper = Zspts(iZ(1)); iZupper = iZ(1);
      Zlower = Zspts(iZ(2)); iZlower = iZ(2);
      if Zupper<Zlower
        Zupper = Zspts(iZ(2)); iZupper = iZ(2);
        Zlower = Zspts(iZ(1)); iZlower = iZ(1);
      end
    end

    % Find the two closest temperature setpoints: "Tupper" and "Tlower"
    dT = abs(T-Tspts); [~,iT] = sort(dT);
    Tupper = Tspts; iTupper = iT;
    Tlower = Tspts; iTlower = iT;
    if length(iT)>1
      Tupper = Tspts(iT(1)); iTupper = iT(1);
      Tlower = Tspts(iT(2)); iTlower = iT(2);
      if Tupper<Tlower
        Tupper = Tspts(iT(2)); iTupper = iT(2);
        Tlower = Tspts(iT(1)); iTlower = iT(1);
      end
    end

    % Compute the four outputs at the time sample
    % Each yk is a column matrix (#TFlocs by 1)
    yk1 = ROMmdls(iTlower,iZlower).C*bigX(:,(iTlower-1)*ZZ+iZlower)...
          + ROMmdls(iTlower,iZlower).D*Iapp;
    yk2 = ROMmdls(iTlower,iZupper).C*bigX(:,(iTlower-1)*ZZ+iZupper)...
          + ROMmdls(iTlower,iZupper).D*Iapp;
    yk3 = ROMmdls(iTupper,iZlower).C*bigX(:,(iTupper-1)*ZZ+iZlower)...
          + ROMmdls(iTupper,iZlower).D*Iapp;
    yk4 = ROMmdls(iTupper,iZupper).C*bigX(:,(iTupper-1)*ZZ+iZupper)...
          + ROMmdls(iTupper,iZupper).D*Iapp;

    % ---------------------------------------------------
    % Step 3: Update states of all pre-computed ROM as
    %         bigX[k+1] = bigA.*bigX[k] + B*iapp[k]
    % ---------------------------------------------------
    bigX = bigA.*bigX + Iapp;
        
    % ----------------------------------------------------------
    % Step 4: Compute the linear output at the present setpoint
    % ----------------------------------------------------------
    % Compute the blending factors
    alphaZ = 0; alphaT = 0;
    if length(iZ)>1
      alphaZ = (ROMout(kks,kkp).cellSOC(k+1)-Zlower)/(Zupper-Zlower); 
    end
    if length(iT)>1
      alphaT = (T-Tlower)/(Tupper-Tlower); 
    end

    % Compute the present linear output yk (a column vector)
    yk = (1-alphaT)*((1-alphaZ)*yk1 + alphaZ*yk2)...
         +alphaT*((1-alphaZ)*yk3 + alphaZ*yk4);

    % -------------------------------------------------------
    % Step 5: Add nonlinear corrections to the linear output
    % -------------------------------------------------------
    % Interfacial total molar rate (ifdl)
    ROMout(kks,kkp).negIfdl(k+1,:) = yk(ROMs.ind.negIfdl);
    ROMout(kks,kkp).posIfdl(k+1,:) = yk(ROMs.ind.posIfdl);
    ROMout(kks,kkp).negIfdl0(k+1)  = yk(ROMs.ind.negIfdl0);
    ROMout(kks,kkp).posIfdl3(k+1)  = yk(ROMs.ind.posIfdl3);

    % Interfacial faradaic molar rate (if)
    ROMout(kks,kkp).negIf(k+1,:)   = yk(ROMs.ind.negIf);
    ROMout(kks,kkp).posIf(k+1,:)   = yk(ROMs.ind.posIf);
    ROMout(kks,kkp).negIf0(k+1)    = yk(ROMs.ind.negIf0);
    ROMout(kks,kkp).posIf3(k+1)    = yk(ROMs.ind.posIf3);

    % Interfacial nonfaradaic molar rate (Idl)
    ROMout(kks,kkp).negIdl(k+1,:)  = yk(ROMs.ind.negIdl);
    ROMout(kks,kkp).posIdl(k+1,:)  = yk(ROMs.ind.posIdl);

    % Solid surface stoichiometries (thetass)
    ROMout(kks,kkp).negThetass(k+1,:) = yk(ROMs.ind.negThetass) + SOC0n;
    if any(ROMout(kks,kkp).negThetass(k+1,:) < 0)
      shortWarn('negThetass < 0'); 
      ROMout(kks,kkp).negThetass(ROMout(kks,kkp).negThetass < 0) = 1e-6;
    end
    if any(ROMout(kks,kkp).negThetass(k+1,:) > 1)
      shortWarn('negThetass > 1'); 
      ROMout(kks,kkp).negThetass(ROMout(kks,kkp).negThetass > 1) = 1-1e-6;
    end

    ROMout(kks,kkp).posThetass(k+1,:) = yk(ROMs.ind.posThetass) + SOC0p;
    if any(ROMout(kks,kkp).posThetass(k+1,:) < 0)
      shortWarn('posThetass < 0'); 
      ROMout(kks,kkp).posThetass(ROMout(kks,kkp).posThetass < 0) = 1e-6;
    end
    if any(ROMout(kks,kkp).posThetass(k+1,:) > 1)
      shortWarn('posThetass > 1'); 
      ROMout(kks,kkp).posThetass(ROMout(kks,kkp).posThetass > 1) = 1-1e-6;
    end

    ROMout(kks,kkp).negThetass0(k+1) = yk(ROMs.ind.negThetass0) + SOC0n;
    if any(ROMout(kks,kkp).negThetass0(k+1) < 0)
      shortWarn('negThetass0 < 0'); 
      ROMout(kks,kkp).negThetass0(ROMout(kks,kkp).negThetass0 < 0) = 1e-6;
    end
    if any(ROMout(kks,kkp).negThetass0(k+1) > 1)
      shortWarn('negThetass0 > 1'); 
      ROMout(kks,kkp).negThetass0(ROMout(kks,kkp).negThetass0 > 1) = 1-1e-6;
    end

    ROMout(kks,kkp).posThetass3(k+1) = yk(ROMs.ind.posThetass3) + SOC0p;
    if any(ROMout(kks,kkp).posThetass3(k+1) < 0)
      shortWarn('posThetass3 < 0'); 
      ROMout(kks,kkp).posThetass3(ROMout(kks,kkp).posThetass3 < 0) = 1e-6;
    end
    if any(ROMout(kks,kkp).posThetass3(k+1) > 1)
      shortWarn('posThetass3 > 1'); 
      ROMout(kks,kkp).posThetass3(ROMout(kks,kkp).posThetass3 > 1) = 1-1e-6;
    end

    % Solid-electrolyte potential difference (phise)
    % The linear output from yk is integrator-removed version
    UocpnAvg = cellData.function.neg.Uocp(ROMout(kks,kkp).negSOC(k+1),T);
    UocppAvg = cellData.function.pos.Uocp(ROMout(kks,kkp).posSOC(k+1),T);
    ROMout(kks,kkp).negPhise(k+1,:) = yk(ROMs.ind.negPhise)  + UocpnAvg;
    ROMout(kks,kkp).posPhise(k+1,:) = yk(ROMs.ind.posPhise)  + UocppAvg;
    ROMout(kks,kkp).negPhise0(k+1)  = yk(ROMs.ind.negPhise0) + UocpnAvg;

    % Compute electrolyte potential: first phie(0,t) then phie(1:3,t)
    ROMout(kks,kkp).Phie(k+1,1)     = 0 - ROMout(kks,kkp).negPhise0(k+1); 
    ROMout(kks,kkp).Phie(k+1,2:end) = yk(ROMs.ind.Phie) - ROMout(kks,kkp).negPhise0(k+1); 

    % Compute electrolyte stoichiometries (thetae)
    ROMout(kks,kkp).Thetae(k+1,:) = yk(ROMs.ind.Thetae) + 1;
    if any(ROMout(kks,kkp).Thetae(k+1,:) < 0)
      shortWarn('Thetae < 0'); ROMout(kks,kkp).Thetae(ROMout(kks,kkp).Thetae < 0) = 1e-6;
    end

    % Compute overpotential at current-collectors via asinh method (eta)
    k0n = cellData.function.neg.k0(ROMout(kks,kkp).negSOC(k+1),T);
    k0p = cellData.function.pos.k0(ROMout(kks,kkp).posSOC(k+1),T);

    i0n = k0n*sqrt(ROMout(kks,kkp).Thetae(k+1,1)* ...
               (1-ROMout(kks,kkp).negThetass0(k+1))*ROMout(kks,kkp).negThetass0(k+1));
    i0p = k0p*sqrt(ROMout(kks,kkp).Thetae(k+1,end)*...
               (1-ROMout(kks,kkp).posThetass3(k+1))*ROMout(kks,kkp).posThetass3(k+1));
    ROMout(kks,kkp).negEta0(k+1) = 2*R*T/F*asinh(ROMout(kks,kkp).negIf0(k+1)/(2*i0n));
    ROMout(kks,kkp).posEta3(k+1) = 2*R*T/F*asinh(ROMout(kks,kkp).posIf3(k+1)/(2*i0p));

    % Compute cell voltage (ROMout(kks,kkp).Vcell)
    Uocpn0 = cellData.function.neg.Uocp(ROMout(kks,kkp).negThetass0(k+1),T);
    Uocpp3 = cellData.function.pos.Uocp(ROMout(kks,kkp).posThetass3(k+1),T);
    Rfn    = cellData.function.neg.Rf(ROMout(kks,kkp).negSOC(k+1),T);
    Rfp    = cellData.function.pos.Rf(ROMout(kks,kkp).posSOC(k+1),T);

    ROMout(kks,kkp).Vcell(k+1) = ROMout(kks,kkp).posEta3(k+1) - ROMout(kks,kkp).negEta0(k+1)...
        + yk(ROMs.ind.Phie(end)) + Uocpp3 - Uocpn0...
        + (Rfp*ROMout(kks,kkp).posIfdl3(k+1) - Rfn*ROMout(kks,kkp).negIfdl0(k+1));

    % Finally, compute solid potential (phis)
    ROMout(kks,kkp).negPhis(k+1,:) = yk(ROMs.ind.negPhis);
    ROMout(kks,kkp).posPhis(k+1,:) = yk(ROMs.ind.posPhis) + ROMout(kks,kkp).Vcell(k+1);
    
    % Update Vcell if Rc is nonzero
    ROMout(kks,kkp).Vcell(k+1) = ROMout(kks,kkp).Vcell(k+1) - Rc*Iapp;
    ROMout(kks,kkp).Iapp(k+1)  = Iapp;
    
    % Function outputs: voltage and updated cell state
    Vcell = ROMout(kks,kkp).Vcell(k+1);
    newCellState.bigX = bigX;
    newCellState.SOCnAvg = SOCnAvg;
    newCellState.SOCpAvg = SOCpAvg;
  end

  %% ----------------------------------------------------------------------
  % This function packages all simulation data into the output data
  % structure to return to the user
  % -----------------------------------------------------------------------
  function storeData
    % A final information update
    fprintf('%d %2.8f \n Voltage: %2.8f\n',length(ik),...
            ROMout(1,1).cellSOC(end),ROMout(1,1).Vcell(end));

    for kks = 1:Ns
      for kkp = 1:Np
        % Save basic information
        ROMout(kks,kkp).blending = 'output-blending';
        ROMout(kks,kkp).cellData = ROMs.ROM(kks,kkp).cellData;
        ROMout(kks,kkp).time     = tk(:);
        % ROMout(kks,kkp).Iapp     = ik(:);
        ROMout(kks,kkp).T        = Tk(:);

        % Save each variable
        ROMout(kks,kkp).Ifdl    = [ROMout(kks,kkp).negIfdl ROMout(kks,kkp).posIfdl];
        ROMout(kks,kkp).If      = [ROMout(kks,kkp).negIf ROMout(kks,kkp).posIf];
        ROMout(kks,kkp).Idl     = [ROMout(kks,kkp).negIdl ROMout(kks,kkp).posIdl];
        ROMout(kks,kkp).Phis    = [ROMout(kks,kkp).negPhis ROMout(kks,kkp).posPhis];
        ROMout(kks,kkp).Phise   = [ROMout(kks,kkp).negPhise ROMout(kks,kkp).posPhise];
        ROMout(kks,kkp).Thetass = [ROMout(kks,kkp).negThetass ROMout(kks,kkp).posThetass];

        % Save locations of each variable
        tfLocs                  = ROMs.tfLocs;
        ROMout(kks,kkp).xLocs.Ifdl    = [tfLocs(ROMs.ind.negIfdl);...
                                tfLocs(ROMs.ind.posIfdl)];
        ROMout(kks,kkp).xLocs.If      = [tfLocs(ROMs.ind.negIf);...
                                tfLocs(ROMs.ind.posIf)];
        ROMout(kks,kkp).xLocs.Idl     = [tfLocs(ROMs.ind.negIdl);...
                                tfLocs(ROMs.ind.posIdl)];
        ROMout(kks,kkp).xLocs.Phis    = [tfLocs(ROMs.ind.negPhis);...
                                tfLocs(ROMs.ind.posPhis)];
        ROMout(kks,kkp).xLocs.Phie    = [0;tfLocs(ROMs.ind.negPhie);...
                                tfLocs(ROMs.ind.sepPhie);...
                                tfLocs(ROMs.ind.posPhie)];
        ROMout(kks,kkp).xLocs.Phise   = [tfLocs(ROMs.ind.negPhise);...
                                tfLocs(ROMs.ind.posPhise)];
        ROMout(kks,kkp).xLocs.Thetass = [tfLocs(ROMs.ind.negThetass);...
                                tfLocs(ROMs.ind.posThetass)];
        ROMout(kks,kkp).xLocs.Thetae  = [tfLocs(ROMs.ind.negThetae);...
                                tfLocs(ROMs.ind.sepThetae);...
                                tfLocs(ROMs.ind.posThetae)];
      end
    end
  end

  %% ----------------------------------------------------------------------
  % This function displays a short warning message to the user if warnState
  % is []. (Does not display line numbers of warning, etc.)
  % -----------------------------------------------------------------------
  function shortWarn(msg)
    if isfield(simData,'warnOff'), return; end
    persistent warnState
    if strcmpi(msg,'on')
      warnState = []; 
    elseif strcmpi(msg,'off')
      warnState = 1;
    elseif isempty(warnState)
      cprintf([1,1/2,0],[' - Warning: ' msg '\n']);
    end
  end
end