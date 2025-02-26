% function [C,Lambda,J,Z,Rct] = tfCommon(s,cellData)
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              transfer functions (TFs) are to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   C        = a matrix of c1n,c2n,c3n,c4n,c1s,c11s,c2s,c1p,c2p,c3p,c4p
%              coefficients required to evaluate TFs 
%   Lambda   = a matrix of the Lambda values required to evaluate TFs
%   J        = a matrix of the "j" values (e.g., required to evaluate ifdl)
%   Z        = a matrix of interface impedances
%   Rct      = a matrix of charge-transfer resistances
%
% This utility function evaluates terms that are common to implementing
% most of the cell transfer functions. It is unlikely that this function
% will ever need to be invoked by the user. It is invoked directly by the
% TF function routines that require its outputs.
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

function [C,Lambda,J,Z,Rct] = tfCommon(s,cellData)
  % first, check to see if we have already computed the relevant data...
  % don't recompute unless necessary!
  if isfield(cellData,'common')
    if isfield(cellData.common,'s')
      if isequal(cellData.common.s,s)
        C = cellData.common.C;
        Lambda = cellData.common.L;
        J = cellData.common.J;
        Z = cellData.common.Z;
        Rct = cellData.common.Rct;
        return
      end
    end
  end

  % Set up negative- and positive-electrode constants
  F = cellData.const.F;
  R = cellData.const.R;
  T = cellData.const.T;
  f = F/(R*T);
  ind0 = find(s==0);

  % Set up constant values
  psi = cellData.const.psi;
  kappaD = cellData.const.kD;
  Q = cellData.const.Q;  

  % set up negative-electrode constants...
  socn = cellData.neg.soc;
  sigman = cellData.neg.sigma;
  kappan = cellData.neg.kappa;
  qen = cellData.neg.qe;
  Dsn = cellData.neg.Ds;
  DeltaQn = abs(cellData.neg.theta100 - cellData.neg.theta0);
  csmaxn = 10800*Q*Dsn/DeltaQn;

  Rf = cellData.neg.Rf; 
  % beta = sqrt(s.^cellData.neg.nF/Dsn);
  beta = sqrt((s/Dsn).^cellData.neg.nF); % 20220914
  Zsn = (cellData.neg.dUocp/csmaxn)...
       *((beta.^2 + 3*(1-beta.*coth(beta)))...
       ./(beta.^2.*(1-beta.*coth(beta))) - 3*Dsn./s);
  Zsn(ind0) = Inf;

  k0 = cellData.neg.k0;
  alpha = cellData.neg.alpha; % charge-transfer alpha
  if length(k0) > 1 % Then, this is an MSMR-kinetics model
    Uref = cellData.neg.Uocp;
    X    = cellData.neg.X;
    U0   = cellData.neg.U0;
    W    = cellData.neg.omega;
    % Compute fractions xj and i0 of each gallary (vectors)
    x    = X./(1+exp(f*(Uref-U0)./W)); % sum(xjn) = thetan
    i0   = sum(k0.*((x).^(W.*alpha)).*((X-x).^(W.*(1-alpha))));  
  else % This is a standard Butler-Volmer-kinetics model
    i0   = k0.*((1 - socn)).^(1 - alpha).*socn.^(alpha);   
  end
  Rctn = R*T/i0/F;
  
  Rdln = cellData.neg.Rdl;
  Cdln = cellData.neg.Cdl;
  nDLn = cellData.neg.nDL;
  wDLn = cellData.neg.wDL;
  Zdln = Rdln + (Cdln*s+wDLn).^(1-nDLn)./(Cdln^(2-nDLn)*s);
  Zsen = Rf + 1./(1./(Rctn + Zsn) + 1./Zdln); 
  Isplitn = 1./(1+(Rctn + Zsn)./Zdln);
  
  % set up separator constants
  kappas = cellData.sep.kappa;
  qes = cellData.sep.qe;

  % set up positive-electrode constants...
  socp = cellData.pos.soc;
  sigmap = cellData.pos.sigma;
  kappap = cellData.pos.kappa;
  qep = cellData.pos.qe;
  Dsp = cellData.pos.Ds;
  DeltaQp = abs(cellData.pos.theta100 - cellData.pos.theta0);
  csmaxp = 10800*Q*Dsp/DeltaQp;

  Rf = cellData.pos.Rf; 
  % beta = sqrt(s.^cellData.pos.nF/Dsp);
  beta = sqrt((s/Dsp).^cellData.pos.nF); % 20220914
  Zsp = (cellData.pos.dUocp/csmaxp)...
       *((beta.^2 + 3*(1-beta.*coth(beta)))...
       ./(beta.^2.*(1-beta.*coth(beta))) - 3*Dsp./s);
  Zsp(ind0) = Inf;

  k0 = cellData.pos.k0;
  alpha = cellData.pos.alpha; % charge-transfer alpha
  if length(k0) > 1 % Then, this is an MSMR-kinetics model
    Uref = cellData.pos.Uocp;
    X    = cellData.pos.X;
    U0   = cellData.pos.U0;
    W    = cellData.pos.omega;
    % Compute fractions xj and i0 of each gallary (vectors)
    x    = X./(1+exp(f*(Uref-U0)./W)); % sum(xjn) = thetan
    i0   = sum(k0.*((x).^(W.*alpha)).*((X-x).^(W.*(1-alpha))));  
  else % This is a standard Butler-Volmer-kinetics model
    i0   = k0.*((1 - socp)).^(1 - alpha).*socp.^(alpha);   
  end
  Rctp = R*T/i0/F;
  
  Rdlp = cellData.pos.Rdl;
  Cdlp = cellData.pos.Cdl;
  nDLp = cellData.pos.nDL;  
  wDLp = cellData.pos.wDL;
  
  Zdlp = Rdlp + (Cdlp*s+wDLp).^(1-nDLp)./(Cdlp^(2-nDLp)*s);
  Zsep = Rf + 1./(1./(Rctp + Zsp) + 1./Zdlp); 
  Isplitp = 1./(1+(Rctp + Zsp)./Zdlp);

  % Set up solution variables
  Lambda1s = sqrt((3600*qes/(psi*T)/kappas)*s); 
  if ~isempty(ind0), Lambda1s(ind0)=0; end

  mu1n = (1/kappan/(psi*T))./Zsen;
  mu2n = (kappaD*T/(psi*T)/kappan)./Zsen - (3600*qen/(psi*T)/kappan).*s;
  tau1n = (1/sigman + 1/kappan)./Zsen - mu2n; 
  tau2n = (3600*qen/(psi*T)/kappan)*s .* (1/sigman + 1/kappan)./Zsen; 
  Lambda1n = sqrt(0.5*(tau1n-sqrt(tau1n.^2-4*tau2n))); 
  Lambda2n = sqrt(0.5*(tau1n+sqrt(tau1n.^2-4*tau2n))); 
  lambda1n = Lambda1n.^3+mu2n.*Lambda1n;
  lambda2n = Lambda2n.^3+mu2n.*Lambda2n; 
  wplusn  = kappas*Lambda1s.*(lambda2n.*coth(Lambda1n)- ...
            lambda1n.*coth(Lambda2n))+kappan*(lambda2n.*Lambda1n- ...
            lambda1n.*Lambda2n); 
  wminusn = kappas*Lambda1s.*(lambda2n.*coth(Lambda1n)- ...
            lambda1n.*coth(Lambda2n))-kappan*(lambda2n.*Lambda1n- ...
            lambda1n.*Lambda2n); 
  if ~isempty(ind0)
    Lambda1n(ind0)=0; Lambda2n(ind0)=0; 
    wplusn(ind0)=0; wminusn(ind0)=0; 
  end

  mu1p = (1/kappap/(psi*T))./Zsep;
  mu2p = (kappaD*T/(psi*T)/kappap)./Zsep - (3600*qep/(psi*T)/kappap).*s;
  tau1p = (1/sigmap + 1/kappap)./Zsep - mu2p; 
  tau2p = (3600*qep/(psi*T)/kappap)*s .* (1/sigmap + 1/kappap)./Zsep; 
  Lambda1p = sqrt(0.5*(tau1p-sqrt(tau1p.^2-4*tau2p))); 
  Lambda2p = sqrt(0.5*(tau1p+sqrt(tau1p.^2-4*tau2p))); 
  lambda1p = Lambda1p.^3+mu2p.*Lambda1p;
  lambda2p = Lambda2p.^3+mu2p.*Lambda2p;
  wplusp  = kappas*Lambda1s.*(lambda2p.*coth(Lambda1p)- ...
            lambda1p.*coth(Lambda2p))+kappap*(lambda2p.*Lambda1p- ...
            lambda1p.*Lambda2p); 
  wminusp = kappas*Lambda1s.*(lambda2p.*coth(Lambda1p)- ...
            lambda1p.*coth(Lambda2p))-kappap*(lambda2p.*Lambda1p- ...
            lambda1p.*Lambda2p); 
  if ~isempty(ind0) 
    Lambda1p(ind0)=0; Lambda2p(ind0)=0; 
    wplusp(ind0)=0; wminusp(ind0)=0; 
  end

  denn = sigman*kappan*(wminusn.*wminusp.*exp(-2*Lambda1s)-wplusn.*wplusp);
  denp = sigmap*kappap*(wminusn.*wminusp.*exp(-2*Lambda1s)-wplusn.*wplusp);

  c1s = -kappan*mu1n.*wminusp.*exp(-Lambda1s)./denn.*(...
        sigman*(Lambda2n.*coth(Lambda1n)-Lambda1n.*coth(Lambda2n)) ...
        +kappan*(Lambda2n.*csch(Lambda1n)-Lambda1n.*csch(Lambda2n))) ...
        +kappap*mu1p.*wplusn./denp.*(...
        sigmap*(Lambda2p.*coth(Lambda1p)-Lambda1p.*coth(Lambda2p)) ...
        +kappap*(Lambda2p.*csch(Lambda1p)-Lambda1p.*csch(Lambda2p)));
  c2s = -kappan*mu1n.*wplusp./denn.*(...
        sigman*(Lambda2n.*coth(Lambda1n)-Lambda1n.*coth(Lambda2n))...
        +kappan*(Lambda2n.*csch(Lambda1n)-Lambda1n.*csch(Lambda2n))) ...
        +kappap*mu1p.*wminusn.*exp(-Lambda1s)./denp.*(...
        sigmap*(Lambda2p.*coth(Lambda1p)-Lambda1p.*coth(Lambda2p))...
        +kappap*(Lambda2p.*csch(Lambda1p)-Lambda1p.*csch(Lambda2p)));
            
  c1n =  1./(kappan*(lambda2n.*Lambda1n-lambda1n.*Lambda2n).* ...
         (1-exp(-2*Lambda1n))).*(mu1n.*Lambda2n.* ...
         (1+kappan*exp(-Lambda1n)/sigman)+kappas*Lambda1s.* ...
         lambda2n.*(c1s.*exp(-Lambda1s)-c2s)); 
  c2n =  1./(kappan*(lambda2n.*Lambda1n-lambda1n.*Lambda2n))...
         .*(mu1n.*Lambda2n.*(csch(Lambda1n)/2+kappan/sigman./ ...
         (1-exp(-2*Lambda1n)))+kappas*Lambda1s.*lambda2n.* ...
         (c1s.*exp(-Lambda1s)-c2s).*csch(Lambda1n)/2);
  c3n = -1./(kappan*(lambda2n.*Lambda1n-lambda1n.*Lambda2n).* ...
        (1-exp(-2*Lambda2n))).*(mu1n.*Lambda1n.* ...
        (1+kappan*exp(-Lambda2n)/sigman)+kappas*Lambda1s.* ...
        lambda1n.*(c1s.*exp(-Lambda1s)-c2s));
  c4n = -1./(kappan*(lambda2n.*Lambda1n-lambda1n.*Lambda2n))...
         .*(mu1n.*Lambda1n.*(csch(Lambda2n)/2+kappan/sigman./ ...
         (1-exp(-2*Lambda2n)))+kappas*Lambda1s.*lambda1n.* ...
         (c1s.*exp(-Lambda1s)-c2s).*csch(Lambda2n)/2);
       
  c1p = -1./(kappap*(lambda2p.*Lambda1p-lambda1p.*Lambda2p).* ...
        (1-exp(-2*Lambda1p))).*(mu1p.*Lambda2p.* ...
        (1+kappap*exp(-Lambda1p)/sigmap)+kappas*Lambda1s.*lambda2p.* ...
        (c1s-c2s.*exp(-Lambda1s)));
  c2p = -1./(kappap*(lambda2p.*Lambda1p-lambda1p.*Lambda2p))...
        .*(mu1p.*Lambda2p.*(csch(Lambda1p)/2+kappap/sigmap./ ...
        (1-exp(-2*Lambda1p)))+kappas*Lambda1s.*lambda2p.* ...
        (c1s-c2s.*exp(-Lambda1s)).*csch(Lambda1p)/2);
  c3p = 1./(kappap*(lambda2p.*Lambda1p-lambda1p.*Lambda2p).* ...
        (1-exp(-2*Lambda2p))).*(mu1p.*Lambda1p.* ...
        (1+kappap*exp(-Lambda2p)/sigmap)+kappas*Lambda1s.*lambda1p.* ...
        (c1s-c2s.*exp(-Lambda1s)));
  c4p = 1./(kappap*(lambda2p.*Lambda1p-lambda1p.*Lambda2p))...
        .*(mu1p.*Lambda1p.*(csch(Lambda2p)/2+kappap/sigmap./ ...
        (1-exp(-2*Lambda2p)))+kappas*Lambda1s.*lambda1p.* ...
        (c1s-c2s.*exp(-Lambda1s)).*csch(Lambda2p)/2);

  j1n = c1n.*(-psi*T*kappan*Lambda1n.^2+3600*qen*s);
  j2n = c2n.*(-psi*T*kappan*Lambda1n.^2+3600*qen*s);
  j3n = c3n.*(-psi*T*kappan*Lambda2n.^2+3600*qen*s);
  j4n = c4n.*(-psi*T*kappan*Lambda2n.^2+3600*qen*s);
  j1p = c1p.*(-psi*T*kappap*Lambda1p.^2+3600*qep*s);
  j2p = c2p.*(-psi*T*kappap*Lambda1p.^2+3600*qep*s);
  j3p = c3p.*(-psi*T*kappap*Lambda2p.^2+3600*qep*s);
  j4p = c4p.*(-psi*T*kappap*Lambda2p.^2+3600*qep*s);

  C = [c1n;c2n;c3n;c4n;c1s;c2s;c1p;c2p;c3p;c4p];    
  Lambda = [Lambda1n;Lambda2n;Lambda1s;Lambda1p;Lambda2p]; 
  J = [j1n;j2n;j3n;j4n;j1p;j2p;j3p;j4p];
  Z = [Zsen; Zsep; Zsn; Zsp; Isplitn; Isplitp]; 
  Rct = [Rctn; Rctp];
end