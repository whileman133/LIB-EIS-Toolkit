% This example shows how to plot the output of a FOM.
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
function plotFOM(FOMout)
  t = FOMout.time/60; % [mins]
  locations = {'neg:cc','neg:sep','pos:sep','pos:cc'};
  M = 1;
  
  % Choose which to plot
  Iapp    = 1;
  Vcell   = 1;

  Phis    = 1;
  Phie    = 1;
  Phise   = 1;
  Thetae  = 1;
  Thetass = 1;
  Ifdl    = 1;
  If      = 1;
  Idl     = 1;

  %% ---------- Input current and output voltage ------------------
  if Iapp == 1
    figure;clf;
    plot(t(1:M:end),FOMout.Iapp(1:M:end));hold on;

    xlabel('Time (min)');ylabel('Current (A)');
    title('Cell current profile');
    grid on; thesisFormat([0.1 0 0.05 0]);
  end

  if Vcell == 1
    figure;clf;
    plot(t(1:M:end),FOMout.Vcell(1:M:end));hold on;

    xlabel('Time (min)');ylabel('Voltage (V)');
    title('Cell terminal voltage');
    grid on; thesisFormat([0.1 0 0.05 0]);
  end

  %% ------------- Find the four locations in FOM data -------------
  ind = zeros(4,1);

  %% ------------------ Phis ------------------
  if Phis == 1
    xLocs = FOMout.xLocs.Phis; 
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:2
      figure;clf;
      plot(t(1:M:end),FOMout.Phis(1:M:end,ind(x+1))); hold on;    
      xlabel('Time (min)');ylabel('Potential (V)');
      title(sprintf('phis:%s',locations{x+1}));
      grid on; thesisFormat([0.2 0 0.05 0]);
    end
  end

  %% ------------------ Phie ------------------
  if Phie == 1
    xLocs = FOMout.xLocs.Phie; 
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:4
      figure;clf;
      plot(t(1:M:end),FOMout.Phie(1:M:end,ind(x))); hold on;

      xlabel('Time (min)');ylabel('Potential (V)');
      title(sprintf('phie:%s',locations{x}));
      grid on; thesisFormat([0.2 0 0.05 0]);
    end
  end

  %% ------------------ Phise ------------------
  if Phise == 1    
    xLocs = FOMout.xLocs.Phise;
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:4
      figure;clf;
      plot(t(1:M:end),FOMout.Phise(1:M:end,ind(x))); hold on;

      xlabel('Time (min)');ylabel('Potential (V)');
      title(sprintf('phise:%s',locations{x}));
      grid on; thesisFormat([0.2 0 0.05 0]);
    end
  end

  %% ------------------ Thetae ------------------
  if Thetae == 1
    xLocs = FOMout.xLocs.Thetae;
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:4
      figure;clf;
      plot(t(1:M:end),FOMout.Thetae(1:M:end,ind(x))); hold on;

      xlabel('Time (min)');ylabel('Concentration ratio (u/l)');
      title(sprintf('thetae:%s',locations{x}));
      grid on; thesisFormat([0.2 0 0.05 0]);
    end
  end

  %% ------------------ Thetass ------------------
  if Thetass == 1
    xLocs = FOMout.xLocs.Thetass;
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:4
      figure;clf;
      plot(t(1:M:end),FOMout.Thetass(1:M:end,ind(x))); hold on;

      xlabel('Time (min)');ylabel('Concentration ratio (u/l)');
      title(sprintf('thetass:%s',locations{x}));
      grid on; thesisFormat([0.2 0 0.05 0]);
    end
  end

  %% ------------------ Ifdl ------------------
  if Ifdl == 1
    xLocs = FOMout.xLocs.Ifdl;
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:4
      figure;clf;
      plot(t(1:M:end),FOMout.Ifdl(1:M:end,ind(x))*1e6); hold on;

      xlabel('Time (min)');ylabel('Total flux (A)');
      title(sprintf('ifdl:%s',locations{x}));
      grid on; thesisFormat([0.2 0 0.05 0]);
    end
  end

  %% ------------------ If ------------------
  if If == 1
    xLocs = FOMout.xLocs.If;
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:4
      figure;clf;
      plot(t(1:M:end),FOMout.If(1:M:end,ind(x))*1e6); hold on;

      xlabel('Time (min)');ylabel('Faradaic flux (A)');
      title(sprintf('if:%s',locations{x}));
      grid on; thesisFormat([0.15 0 0.05 0]);
    end
  end

  %% ------------------ Idl ------------------
  if Idl == 1
    xLocs = FOMout.xLocs.Idl;
    ind(1) = find(xLocs == 0,1); % negative current-collector
    ind(2) = find(xLocs == 1,1); % neg/sep boundary
    ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
    ind(4) = find(xLocs == 3,1); % positive current-collector

    for x = 1:4
      figure;clf;
      plot(t(1:M:end),FOMout.Idl(1:M:end,ind(x))*1e6); hold on;

      xlabel('Time (min)');ylabel('Double-layer flux (A)');
      title(sprintf('idl:%s',locations{x}));
      grid on; thesisFormat([0.15 0 0.05 0]);
    end
  end
end        