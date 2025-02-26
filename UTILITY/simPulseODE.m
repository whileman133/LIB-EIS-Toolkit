% ---------------------------------------------------------------------------------
% This code applys ODE solver (bvp5c) to simulate pulse response of the cell
%
% Assumption: when the input is a current pulse (short duration and high mag.),
%             concentration variables cannot change instantly. Only potential 
%             variables and interface flux terms can change immediately.
% Mathematical model backgrounds refer to the course ECE5730, Section 3.3
%
% Simulate ODE for phise. Let  x = [0,1]; y = [phisen;phisen';phisep;phisep'] 
%
% Function input:   socvals   (soc stepoints)
%                   cratevals (input current Crate)
%                   p         (cell parameters)
% Function output:  R0 (cell impedance under all soc setpoints and Crates)
% ---------------------------------------------------------------------------------

% THIS IS THE FASTER VERSION THAT HANDLES ALL FOUR CELL LOCS AT ONCE

function R0 = simPulseODE(socvals,cratevals,p)
  %-----------------------------------
  % Unpack p into variables
  %-----------------------------------
  % define fields in "p" that we must keep; delete all others
%   keepFields = {'F';'R';'T';'Q';'alphan';'alphap';'k0n';'k0p';'sigman';...
%       'sigmap';'kappan';'kappas';'kappap';'Rdln';'Rdlp';'Rfn';'Rfp';...
%       'theta0n';'theta0p';'theta100n';'theta100p'};
  keepFields = {'F';'R';'T';'Q';'alphan';'alphap';'k0n';'k0p';...
    'sigman';'sigmap';'kappan';'kappas';'kappap';...
    'Rdln';'Rdlp';'Rfn';'Rfp';...
    'theta0n';'theta0p';'theta100n';'theta100p';...
    'Uocpn';'Uocpp';'U0n';'Xn';'Wn';'U0p';'Xp';'Wp'};    
  p = rmfield(p,setdiff(fieldnames(p),keepFields));

  % ensure that we actually have all the fields that we need
  missingFields = setdiff(keepFields,fieldnames(p));
  if ~isempty(missingFields)
      fprintf('Missing field in parameter vector p: %s\n', missingFields{:});
      error('Cannot continue with missing fields.');
  end

  % sort fields in alphabetic order, "deal" to sorted list of variables
  P = struct2cell(orderfields(p));
%   [F,Q,R,Rdln,Rdlp,Rfn,Rfp,T,alphan,alphap,k0n,k0p,kappan,kappap,kappas,...
%       sigman,sigmap,theta0n,theta0p,theta100n,theta100p] = deal(P{:});
  [F,Q,R,Rdln,Rdlp,Rfn,Rfp,T,U0n,U0p,Uocpn,Uocpp,Wn,Wp,Xn,Xp,...
    alphan,alphap,k0n,k0p,kappan,kappap,kappas,...
    sigman,sigmap,theta0n,theta0p,theta100n,theta100p] = deal(P{:});

  %------------------------------------------------------------
  % Compute R0 using bvp5c solver (solve for phise equation)
  %------------------------------------------------------------
  % Create storage
  dv   = zeros(length(socvals),length(cratevals));
  iapp = zeros(length(socvals),length(cratevals));

  if Rfn == 0 && Rdln == 0
    warning('Rdl and Rf (neg) cannot both be zero');
    R0 = NaN(size(dv));
    return
  end
  if Rfp == 0 && Rdlp == 0
    warning('Rdl and Rf (pos) cannot both be zero');
    R0 = NaN(size(dv));
    return
  end

  % A little faster with "Vectorize" turned on
  bvpOpts = bvpset('Stats','off','Vectorized','on'); % options

  f   = F/(R*T);
  K1n = (1-alphan)*f;
  K2n = -alphan*f;
  K1p = (1-alphap)*f;
  K2p = -alphap*f;

  for kz = 1:length(socvals)
    z = socvals(kz); % the SOC setpoint
    thetan = theta0n + z*(theta100n-theta0n);
    thetap = theta0p + z*(theta100p-theta0p);
    
    if length(k0n) > 1 % Then, this is an MSMR-kinetics model
      Urefn = Uocpn(thetan,T); 
      Urefp = Uocpp(thetap,T); 

      % Compute fractions xj and i0 of each gallary (vectors)
      xn = Xn./(1+exp(f*(Urefn-U0n)./Wn)); % sum(xjn) = thetan
      xp = Xp./(1+exp(f*(Urefp-U0p)./Wp)); % sum(xjp) = thetap
      i0n = k0n.*((xn).^(Wn.*alphan)).*((Xn-xn).^(Wn.*(1-alphan)));
      i0p = k0p.*((xp).^(Wp.*alphap)).*((Xp-xp).^(Wp.*(1-alphap)));
    else % This is a standard Butler-Volmer-kinetics model
      i0n = k0n.*((1-thetan)^(1-alphan))*(thetan^alphan);
      i0p = k0p.*((1-thetap)^(1-alphap))*(thetap^alphap);
    end
    negFn = @(x,phise) sum(i0n.*(exp(K1n.*(phise-Rfn*x))... % x = ifdln
        -exp(K2n.*(phise-Rfn*x)))) + (phise-Rfn*x)/Rdln - x;
    posFn = @(x,phise) sum(i0p.*(exp(K1p.*(phise-Rfp*x))... % x = ifdlp
        -exp(K2p.*(phise-Rfp*x)))) + (phise-Rfp*x)/Rdlp - x;

    negMean = @(x,in) x/Rdln+sum(i0n.*(exp(K1n.*x)-exp(K2n.*x)))-in; % x = etan
    posMean = @(x,ip) x/Rdlp+sum(i0p.*(exp(K1p.*x)-exp(K2p.*x)))+ip; % x = etap
    
    for kc = 1:length(cratevals)
      ik = cratevals(kc)*Q; % applied current (A)
      iapp(kz,kc) = ik;

      if Rdln == 0
        phisemean = Rfn*ik;
      else
        phisemean = fzeroFaster(negMean,0,[],ik) + Rfn*ik;
      end

      if Rdlp == 0
        phisemeap = Rfp*ik;
      else
        phisemeap = fzeroFaster(posMean,0,[],ik) + Rfp*ik;
      end

      % Define initial values [phisen;phisen';phisep;phisep']
      bn = -ik/sigman; an = (+ik/kappan-bn)/2; cn = phisemean-an/3-bn/2;
      bp = +ik/sigmap; ap = (-ik/kappap-bp)/2; cp = phisemeap-ap/3-bp/2;
      yinit = @(x) [an*x.^2+bn*x+cn;2*an*x+bn;...
                    ap*(1-x).^2+bp*(1-x)+cp; ap*(2*(1-x)-6)-bp];
      init.solver = 'bvpinit';
      init.x      = [0 0.5 1];
      init.y      = yinit([0 0.5 1]);
      init.yinit  = yinit;

      sol5c = bvp5c(@odefun,@bcfun,init,bvpOpts); % bvp5c ODE solver
      phisen = [sol5c.y(1,1) sol5c.y(1,end)]; % phisen(0) and phisen(1)
      phisep = [sol5c.y(3,1) sol5c.y(3,end)]; % phisep(2) and phisep(3)

      % ------------------------------------------------------
      %  Compute dv
      % ------------------------------------------------------
      dv(kz,kc) = -ik/(kappan+sigman) -ik/kappas -ik/(kappap+sigmap) ...
          - sigman/(kappan+sigman)*phisen(1,2) ... % phisen(xtilde=1)
          - kappan/(kappan+sigman)*phisen(1,1) ... % phisen(xtilde=0)
          + kappap/(kappap+sigmap)*phisep(1,2) ... % phisep(xtilde=3)
          + sigmap/(kappap+sigmap)*phisep(1,1);    % phisep(xtilde=2)
    end
  end

  % Compute R0
  R0 = -dv./iapp;

  % ----------------------------------------------------------------------
  % Define ODE function: "x" input is not used
  % y    = [phisen;phisen';phisep;phisep'] 
  % dydx = [phisen';phisen'';phisep';phisep'']
  % ----------------------------------------------------------------------
  function dydx = odefun(~,y)
    % Define dydx and compute phisen' and phisep'
    dydx = 0*y; dydx(1,:) = y(2,:); dydx(3,:) = y(4,:);

    % Compute for phisen'' first
    if Rfn == 0 && Rdln ~= 0 % faster solution: 9/4/2020
      dydx(2,:) = (1/sigman+1/kappan)...
          *(i0n*(exp(K1n*y(1,:))-exp(K2n*y(1,:)))...
          +y(1,:)/Rdln);
    elseif Rfn == 0 && Rdln == 0
      warning('Rdl and Rf (neg) cannot both be zero');
      dydx(2,:) = NaN;
    elseif Rfn ~=0 && Rdln == 0
      dydx(2,:) = (1/sigman+1/kappan)*y(1,:)/Rfn;
    else % Rfn ~= 0 && Rdln ~= 0
      guessN = ik;
      for k = 1:size(y,2)
        % Solve "ifdl" from kinetics equation: ifdl = if + idl
        ifdl = fzeroFaster(negFn,guessN,[],y(1,k));
        dydx(2,k) = (1/sigman+1/kappan)*ifdl;
        guessN = ifdl;
      end
      % Tried "fsolve" to vectorize. It worked, but was slower
      % guessN = i*ones(1,size(y,2));
      % ifdl = fsolve(negFn,guessN,[],y(1,:));
      % dydx(2,:) = (1/sigman + 1/kappan)*ifdl;
    end

    % Compute for phisep'' next
    if Rfp == 0 && Rdlp ~= 0 % faster solution: 9/4/2020
      dydx(4,:) = (1/sigmap+1/kappap)...
          *(i0p*(exp(K1p*y(3,:))-exp(K2p*y(3,:)))...
          +y(3,:)/Rdlp);
    elseif Rfp == 0 && Rdlp == 0
      warning('Rdl and Rf (pos) cannot both be zero');
      dydx(4,:) = NaN;
    elseif Rfp ~= 0 && Rdlp == 0
      dydx(4,:) = (1/sigmap+1/kappap)*y(3,:)/Rfp;
    else
      guessP = -ik;
      for k = 1:size(y,2)
        ifdl = fzeroFaster(posFn,guessP,[],y(3,k));
        dydx(4,k) = (1/sigmap+1/kappap)*ifdl;
        guessP = ifdl;
      end
    end
  end

  % ------------------------------------------------------------------------
  % Set up boundary conditions
  % y = [phisen;phisen';phisep;phisep'];
  % yl is the left-side boundary, yr is the right-side boundary
  % ------------------------------------------------------------------------
  function res = bcfun(yl,yr)
    res = [yl(2)+ik/sigman;  % xtilde = 0
           yl(4)-ik/kappap;  % xtilde = 2
           yr(2)-ik/kappan;  % xtilde = 1
           yr(4)+ik/sigmap]; % xtilde = 3
  end
end

% -------------------------------------------------------------------------
% Replace MATLAB's built-in FZERO command with this stripped-down version
% that removes a lot of error checking and tracing capabilities in
% preference for speed. Runs about twice as fast. See help on FZERO
% -------------------------------------------------------------------------
function b = fzeroFaster(FunFcn,x,options,varargin) %#ok<INUSL>
  % I always call with four inputs: anon function, initial guess, empty
  % options [], and a second argument to the anonymous function.
  % I always call with one output, so other outputs from FZERO removed.

  % Initialization
  fcount = 0;
  iter = 0;
  intervaliter = 0;
  tol = eps;

  % Put first feval in try catch
  try
    fx = FunFcn(x,varargin{:});
  catch ME
    if ~isempty(Ffcnstr)
      error('MATLAB:fzero:InvalidFunctionSupplied',...
          getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',...
          sprintf('%s ==> %s','function_handle',Ffcnstr),ME.message)));
    else
      error('MATLAB:fzero:InvalidFunctionSupplied',...
          getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',...
          'function_handle',ME.message)));
    end
  end
  fcount = fcount + 1;
  if fx == 0
    b = x;
    return
  elseif ~isfinite(fx) || ~isreal(fx)
    error('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite',...
        getString(message('MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite')));
  end

  if x ~= 0
    dx = x/50;
  else
    dx = 1/50;
  end

  % Find changes of sign.
  twosqrt = sqrt(2);
  a = x; fa = fx; b = x; fb = fx;

  while (fa > 0) == (fb > 0)
    intervaliter = intervaliter + 1;
    dx = twosqrt*dx;
    a = x - dx;  fa = FunFcn(a,varargin{:});
    fcount = fcount + 1;
    if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
        b = NaN; return
    end

    if (fa > 0) ~= (fb > 0) % check for different sign
      break
    end

    b = x + dx;  fb = FunFcn(b,varargin{:});
    if ~isfinite(fb) || ~isreal(fb) || ~isfinite(b)
      b = NaN; return
    end
    fcount = fcount + 1;
  end % while

  fc = fb;
  while fb ~= 0 && a ~= b
    % Ensure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
      c = a;  fc = fa;
      d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
      a = b;    b = c;    c = a;
      fa = fb;  fb = fc;  fc = fa;
    end

    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0)
      break
    end

    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
      % Bisection
      d = m;  e = m;
    else
      % Interpolation
      s = fb/fa;
      if (a == c)
        % Linear interpolation
        p = 2.0*m*s;
        q = 1.0 - s;
      else
        % Inverse quadratic interpolation
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
        q = (q - 1.0)*(r - 1.0)*(s - 1.0);
      end
      if p > 0
        q = -q;
      else
        p = -p;
      end
      % Is interpolated point acceptable
      if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
        e = d;  d = p/q;
      else
        d = m;  e = m;
      end
    end % Interpolation

    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler
      b = b + d;
    elseif b > c
      b = b - toler;
    else
      b = b + toler;
    end
    fb = FunFcn(b,varargin{:});
    fcount = fcount + 1;
    iter = iter + 1;
  end % Main loop
end