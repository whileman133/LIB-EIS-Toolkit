% function cellParams = loadCellParams(fileName)
% 
% Input:
%   fileName   = full path to Excel spreadsheet of cell parameter values
% Output:
%   cellParams = data structure containing all parameter values
%
% This utility function creates a cell-model data structure by loading the
% values from an Excel spreadsheet. See the "Instructions" tab on any of
% the example spreadsheets in the toolbox "CELL_DEFNS" folder for the
% required file format. We recommend that when you are creating the
% definition for a new cell that you begin by duplicating one of these
% example spreadsheets, then overwriting the values in the duplicate
% spreadsheet that must be changed to describe the new cell.
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

% -- Changelog --
% 04.28.2022 | Wesley Hileman | Edited
%   - readParamTable: use `readcell` instead of `xlsread` (fixes sheet 
%     detection issue when reading model created by saveCellParams)
%   - loadTable: use `readmatrix` instead of `xlsread` (fixes sheet 
%     detection issue when reading model created by saveCellParams)
%   - Extract part of readParamTable into externally-accessible 
%     loadParamTable function.

function cellParams = loadCellParams(fileName)
  inputCellData = readParamTable(fileName,'Parameters'); % read Excel sheet
  if ~inputCellData.lumped % standard parameters, convert to lumped
    inputCellData = convertStandardToLumped(inputCellData);
  end
  cellParams = cellStr2Fun(inputCellData,1); % convert strings to functions
  if cellParams.MSMR
    cellParams.function.neg.Uocp  = makeOCP(cellParams.function.neg);
    cellParams.function.pos.Uocp  = makeOCP(cellParams.function.pos);
    cellParams.function.neg.dUocp = makedOCP(cellParams.function.neg);
    cellParams.function.pos.dUocp = makedOCP(cellParams.function.pos);
  end
  cellParams = makeOCV(cellParams);
end

% Read the Excel spreadsheet from file=fileName, tab=sheet
function cellParams = readParamTable(fileName,sheet)
  cellParams = []; % initialize to empty structure
  
  % Read contents of Worksheet into cell array.
  [ind,data] = loadParamTable(fileName,sheet);
  
  % Individually process the main sections of the XLSX data file
  cellParams.const.F = 96485.3365; % Faraday's constant
  cellParams.const.R = 8.3144621;  % Universal gas constant
  cellParams = processGen(data(ind.gen,:),cellParams);
  cellParams = processCst(data(ind.const,:),cellParams,fileName);
  cellParams = processNeg(data(ind.neg,:),cellParams,fileName);
  cellParams = processSep(data(ind.sep,:),cellParams,fileName);
  cellParams = processPos(data(ind.pos,:),cellParams,fileName);

  % Add SOC function to each electrode
  theta0_neg = str2func(cellParams.function.neg.theta0); 
  theta100_neg = str2func(cellParams.function.neg.theta100);
  theta0_pos = str2func(cellParams.function.pos.theta0); 
  theta100_pos = str2func(cellParams.function.pos.theta100);
  cellParams.function.neg.soc = sprintf('@(x,T) (%g + x*(%g))',...
    theta0_neg(),theta100_neg()-theta0_neg());
  cellParams.function.pos.soc = sprintf('@(x,T) (%g + x*(%g))',...
    theta0_pos(),theta100_pos()-theta0_pos());  
end

% Process the "#general" segment of the data file
% ... cell name and type (lumped or standard)
function cellParams = processGen(data,cellParams)
  % Prepare data by extracting only relevant part, deleting "Code Name"
  vars = data(:,2); vals = data(:,3);
  indVar = cellfun(@isstr,vars); % index of strings in "vars" column
  indVar2 = logical(indVar.*indVar([end,1:end-1])); % Delete "Code Name"
  vars = vars(indVar2); vals = vals(indVar2); 
  indName = strcmpi(vars,'name'); % Search for cell-name row
  if sum(indName) == 0
    error('Missing "name" entry in the general section of xlsx file');
  end
  cellParams.name = vals{indName};
  indLumped = strcmpi(vars,'lumped'); % Search for cell-type row
  if isempty(indLumped)
    error('Missing "lumped" entry in the general section of xlsx file');
  end
  cellParams.lumped = vals{indLumped};
  indMSMR = strcmpi(vars,'MSMR'); % Search for kinetics-type row
  if sum(indMSMR) == 0
    cellParams.MSMR = 0;
  else
    cellParams.MSMR = vals{indMSMR};
  end
end

% Process the "#const" segment of the data file
function cellParams = processCst(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    const.(vars{j}) = vals{j};
    if Eact{j} ~= 0 && ~isnan(Eact{j})
      const.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.const = const;
%   cellParams.function.const = addTemps(const);
end

% Process the "#neg" segment of the data file
function cellParams = processNeg(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    neg.(vars{j}) = vals{j};
    if ischar(Eact{j}) % then probably a vector
      EactVect = str2num(Eact{j}); %#ok<ST2NM>
      EactVect = 1000*EactVect(:);
      EactStr = sprintf(' %g;',EactVect);
      neg.(sprintf('%s__Ea',vars{j})) = sprintf('[%s]',EactStr(1:end-1)); % delete ";" from end
    elseif Eact{j} ~= 0 && ~isnan(Eact{j})
      neg.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.neg = neg;
end

% Process the "#sep" segment of the data file
function cellParams = processSep(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    sep.(vars{j}) = vals{j};
    if ischar(Eact{j}) % then probably a vector
      EactVect = str2num(Eact{j}); %#ok<ST2NM>
      EactVect = 1000*EactVect(:);
      EactStr = sprintf(' %g;',EactVect);
      sep.(sprintf('%s__Ea',vars{j})) = sprintf('[%s]',EactStr(1:end-1)); % delete ";" from end
    elseif Eact{j} ~= 0 && ~isnan(Eact{j})
      sep.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.sep = sep;
end

% Process the "#pos" segment of the data file
function cellParams = processPos(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    pos.(vars{j}) = vals{j};
    if ischar(Eact{j}) % then probably a vector
      EactVect = str2num(Eact{j}); %#ok<ST2NM>
      EactVect = 1000*EactVect(:);
      EactStr = sprintf(' %g;',EactVect);
      pos.(sprintf('%s__Ea',vars{j})) = sprintf('[%s]',EactStr(1:end-1)); % delete ";" from end
    elseif Eact{j} ~= 0 && ~isnan(Eact{j})
      pos.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.pos = pos;
end

function [vars,vals,Eact] = parseData(data,fileName)
  % Prepare data by extracting only relevant part, deleting "Code Name"
  vars = data(:,2); vals = data(:,3); Eact = data(:,4);
  indVar = cellfun(@isstr,vars); % Index of strings
  indVar2 = logical(indVar.*indVar([end,1:end-1])); % Delete "Code Name"
  vars = vars(indVar2); vals = vals(indVar2); Eact = Eact(indVar2);
    
  indNum = cellfun(@isnumeric,vals);
  valsNum = sprintf('''@(x) (%g)''\n',[vals{indNum}]'); 
  vals(indNum) = eval(['{',valsNum,'}']);
  % Replace Inf, if any
  % strInf = strfind(vals,'Inf'); indInf = ~cellfun(@isempty,strInf);
  indInf = contains(vals,'Inf');
  vals(indInf) = {'@(x) (Inf)'};
  
  % Load vectors, if any
  vectInd = find(contains(vals,'['));
  if ~isempty(vectInd)
    for k = 1:length(vectInd)
      vectCell = vals(vectInd(k));
      vectCell = str2num(vectCell{:}); %#ok<ST2NM>
      vectCell = vectCell(:); % force to be column vector
      vectStr = sprintf(' %g;',vectCell);
      vals(vectInd(k)) = {['@(x) ([',vectStr(1:end-1),'])']};
    end
  end
  
  % Load tables, if any
  tableInd = find(contains(vals,'#'));
  if ~isempty(tableInd)
    for k = 1:length(tableInd)
      tabCell = vals(tableInd(k));
      tabCell = tabCell{1};
      tabName = tabCell(2:end);
      tableString = loadTable(fileName,tabName);
      vals(tableInd(k)) = {tableString};
    end
  end

  % Define all input arguments to be x and T
  for nf = 1:length(vals)
    fun_str = vals{nf};
    strFun_del1 = strfind(fun_str,'@(');
    strFun_del2 = strfind(fun_str,')');
    ind_del = find(strFun_del2>strFun_del1,1,'first');
    ind_arg = strFun_del1:strFun_del2(ind_del);
    arg_str = fun_str(ind_arg);
    vals{nf} = strrep(fun_str,arg_str,'@(x,T)');   
  end
end

% Combine Arrhenius relationship with reference parameters in function
function reg = addTemps(reg)
  R = 8.3144621;  % Universal gas constant
  f = fields(reg);
  for k = 1:length(f)
    fk = f{k};
    if length(fk)>4 && strcmp(fk(end-3:end),'__Ea')
      base = fk(1:end-4);
      if isfield(reg,base) % okay, we have a <param> and <param>__Ea combo
        baseExt = sprintf('.*exp(%s*(1/298.15-1/T)/%g)',reg.(fk),R);
        reg.(base) = [reg.(base) baseExt];
      end
      reg = rmfield(reg,fk);
    end
  end
end

% Converts from standard parameters to lumped parameters
% Output "cellData" is structure with STRINGS, not FUNCTIONS at this point
function [ cellData ] = convertStandardToLumped(cellDataStd)
  cellData = cellDataStd;
  cellData.standard = cellData.function;
  cellData = rmfield(cellData,'function'); % remove standard functions, replace
  cellDataFunc = cellStr2Fun(cellDataStd,0); % convert strings to functions
  cellDataFunc = cellDataFunc.function;

  % Regions length
  Lneg = cellDataFunc.neg.L();
  Lsep = cellDataFunc.sep.L();
  Lpos = cellDataFunc.pos.L(); 
  % Set up constants and calculations needed by transfer function
  F = cellData.const.F;
  R = cellData.const.R;
  t0plus  = cellDataFunc.const.t0plus();
  RsNeg = cellDataFunc.neg.Rs(); 
  RsPos = cellDataFunc.pos.Rs();
  if isfield(cellDataFunc.neg,'Dsref')
    DsrefNeg = cellDataFunc.neg.Dsref();  
  else
    DsNeg = cellDataFunc.neg.Ds();
  end
  if isfield(cellDataFunc.pos,'Dsref')
    DsrefPos = cellDataFunc.pos.Dsref(); 
  else
    DsPos = cellDataFunc.pos.Ds();
  end
  sEpsNeg  = cellDataFunc.neg.sEps();
  sEpsPos  = cellDataFunc.pos.sEps();
  Acell = cellDataFunc.const.A();
  Rc    = cellDataFunc.const.Rc();
  asNeg = 3*sEpsNeg/RsNeg;
  asPos = 3*sEpsPos/RsPos;
  eEpsNeg  = cellDataFunc.neg.eEps();
  eEpsSep  = cellDataFunc.sep.eEps();
  eEpsPos  = cellDataFunc.pos.eEps();
  De = cellDataFunc.const.De();
  % We could create "effective" De, but because we assume brugKappa and
  % brugDe are the same, the effective De is unused. Instead, we define psi
  % DeNeg  = cellDataFunc.const.De() * eEpsNeg^cellDataFunc.neg.brugDe();
  % DeNep  = cellDataFunc.const.De() * eEpsSep^cellDataFunc.sep.brugDe();
  % DePos  = cellDataFunc.const.De() * eEpsPos^cellDataFunc.pos.brugDe();
  ce0         = cellDataFunc.const.ce0();
  kappa       = cellDataFunc.const.kappa(ce0,0);
  kappaEffNeg = kappa*(eEpsNeg)^(cellDataFunc.neg.brugDeKappa());
  kappaEffSep = kappa*(eEpsSep)^(cellDataFunc.sep.brugDeKappa());
  kappaEffPos = kappa*(eEpsPos)^(cellDataFunc.pos.brugDeKappa());
  sigmaNeg    = cellDataFunc.neg.sigma();
  sigmaEffNeg = sigmaNeg*(sEpsNeg).^cellDataFunc.neg.brugSigma();
  sigmaPos    = cellDataFunc.pos.sigma();
  sigmaEffPos = sigmaPos.*(sEpsPos).^cellDataFunc.pos.brugSigma();
  theta0Neg   = cellDataFunc.neg.theta0();
  theta100Neg = cellDataFunc.neg.theta100();
  theta0Pos   = cellDataFunc.pos.theta0();
  theta100Pos = cellDataFunc.pos.theta100();
  csmaxNeg    = cellDataFunc.neg.csmax();
  alphaNeg    = cellDataFunc.neg.alpha();
  alphaPos    = cellDataFunc.pos.alpha();
  knormNeg    = cellDataFunc.neg.knorm();
  knormPos    = cellDataFunc.pos.knorm();
  RfNeg       = cellDataFunc.neg.Rf();
  RfPos       = cellDataFunc.pos.Rf();
  dlnfdlnc    = cellDataFunc.const.dlnfdlnc();
  RdlNeg      = cellDataFunc.neg.Rdl();
  RdlPos      = cellDataFunc.pos.Rdl();
  nFNeg       = cellDataFunc.neg.nF();
  nFPos       = cellDataFunc.pos.nF();
  nDLNeg      = cellDataFunc.neg.nDL();
  nDLPos      = cellDataFunc.pos.nDL();
  CdlNeg      = cellDataFunc.neg.Cdl();
  CdlPos      = cellDataFunc.pos.Cdl();
  wDLNeg      = cellDataFunc.neg.wDL();
  wDLPos      = cellDataFunc.pos.wDL();

  % Generate Lumped parameters
  package = @(x) sprintf('@(x,T)(%g)',x);
  fn.neg.Uocp      = cellData.standard.neg.Uocp;
  fn.neg.dUocp     = cellData.standard.neg.dUocp;
  fn.neg.sigma     = package(sigmaEffNeg*Acell/Lneg);
  fn.neg.kappa     = package(kappaEffNeg*Acell/Lneg);
  fn.neg.theta0    = package(theta0Neg);
  fn.neg.theta100  = package(theta100Neg);
  if isfield(cellDataFunc.neg,'Dsref') 
    fn.neg.Dsref   = package(DsrefNeg/RsNeg^2);
  else
    fn.neg.Ds      = package(DsNeg/RsNeg^2);
  end
  fn.neg.qe        = package(eEpsNeg*ce0*Acell*Lneg*F/(1-t0plus)/3600);
  fn.neg.k0        = package(knormNeg*asNeg*Acell*Lneg*F);
  fn.neg.alpha     = package(alphaNeg);
  fn.neg.Rf        = package(RfNeg/(asNeg*Acell*Lneg));
  fn.neg.nF        = package(nFNeg);
  fn.neg.Cdl       = package(CdlNeg*asNeg*Acell*Lneg);
  fn.neg.nDL       = package(nDLNeg);
  fn.neg.wDL       = package(wDLNeg);
  fn.neg.Rdl       = package(RdlNeg/(asNeg*Acell*Lneg));
  fn.sep.kappa     = package(kappaEffSep*Acell/Lsep);
  fn.sep.qe        = package(eEpsSep*ce0*Acell*Lsep*F/(1-t0plus)/3600);

  fn.pos.Uocp      = cellData.standard.pos.Uocp;
  fn.pos.dUocp     = cellData.standard.pos.dUocp;
  fn.pos.sigma     = package(sigmaEffPos*Acell/Lpos);
  fn.pos.kappa     = package(kappaEffPos*Acell/Lpos);
  fn.pos.theta0    = package(theta0Pos);
  fn.pos.theta100  = package(theta100Pos);
  if isfield(cellDataFunc.pos,'Dsref') 
    fn.pos.Dsref   = package(DsrefPos/RsPos^2);
  else
    fn.pos.Ds      = package(DsPos/RsPos^2);
  end
  fn.pos.qe        = package(eEpsPos*ce0*Acell*Lpos*F/(1-t0plus)/3600);
  fn.pos.k0        = package(knormPos*asPos*Acell*Lpos*F);
  fn.pos.alpha     = package(alphaPos);
  fn.pos.Rf        = package(RfPos/(asPos*Acell*Lpos));
  fn.pos.nF        = package(nFPos);
  fn.pos.Cdl       = package(CdlPos*asPos*Acell*Lpos);
  fn.pos.nDL       = package(nDLPos);
  fn.pos.wDL       = package(wDLPos);
  fn.pos.Rdl       = package(RdlPos/(asPos*Acell*Lpos));
  fn.const.Q       = package(sEpsNeg*Acell*Lneg*csmaxNeg*...
                           abs(theta100Neg - theta0Neg)*F/3600);
  fn.const.Rc      = package(Rc);
  fn.const.psi     = sprintf('@(x,T)(%g/T)',F*De/kappa*ce0/(1-t0plus));
  fn.const.kD      = package(2*R*(t0plus-1)/F*(1+dlnfdlnc));
  
  % Copy activation energies...
  reg = {'const','neg','sep','pos'};
  for nr = 1:length(reg)
    regData = cellData.standard.(reg{nr});
    f = fields(regData);
    for k = 1:length(f)
      fk = f{k};
      if strcmp(fk,'kappa__Ea') % copy from "const" to other regions
        fn.neg.kappa__Ea = regData.kappa__Ea;
        fn.sep.kappa__Ea = regData.kappa__Ea;
        fn.pos.kappa__Ea = regData.kappa__Ea;
      end
      if length(fk)>4 && strcmp(fk(end-3:end),'__Ea')
        fnReg = fn.(reg{nr});
        if strcmp(fk(1:end-4),'knorm')
          fnReg.k0__Ea = regData.(fk);
        else
          fnReg.(fk) = regData.(fk);
        end
        fn.(reg{nr}) = fnReg;
      end
    end
  end
  
  % copy electrode SOC equation
  fn.neg.soc = cellData.standard.neg.soc;
  fn.pos.soc = cellData.standard.pos.soc;
  cellData.function = fn;
end

function tableStr = loadTable(fileName,tabName)
  tabData = readtable(fileName,'Sheet',tabName,'ReadVariableNames',false);
  tabData = table2array(tabData);
  tabData1 = mat2str(tabData(:,1));  % independent variable
  tabData2 = mat2str(tabData(:,2));  % dependent variable
  % Replace data by string function for table lookup
  tableStr = sprintf('@(x,T) (interp1(%s,%s,x,''linear'',''extrap''))',...
             tabData1,tabData2);
end

% Convert strings stored in the data structure to MATLAB inline functions
function cellData = cellStr2Fun(cellData,temps)
  reg = {'const','neg','sep','pos'};
  for nr = 1:length(reg)
    if temps
      cellData.function.(reg{nr}) = addTemps(cellData.function.(reg{nr}));
    end
    y1 = cellData.function.(reg{nr});
    f1 = fields(y1);
    y2 = struct2cell(y1);
    fnInd = find(contains(y2,'@'));
    y3 = cellfun(@str2func,y2(fnInd),'UniformOutput',false);
    y4 = cell2struct(y3,f1(fnInd));
    cellData.function.(reg{nr}) = y4;
  end
end

% Make cell-level OCV function from electrode-level OCP functions
function cellParams = makeOCV(cellParams)
  SOC = 0:0.001:1;
  negTheta = cellParams.function.neg.soc(SOC,298.15);
  posTheta = cellParams.function.pos.soc(SOC,298.15);
  negOCP = cellParams.function.neg.Uocp(negTheta,298.15);
  posOCP = cellParams.function.pos.Uocp(posTheta,298.15);
  cellOCV = posOCP - negOCP;
  SOCstr = sprintf(' %g;',SOC);
  SOCstr = SOCstr(1:end-1); % remove trailing ';'
  OCVstr = sprintf(' %g;',cellOCV);
  OCVstr = OCVstr(1:end-1); % ditto
  OCVfn = sprintf('@(x,T) interp1([%s],[%s],x,''linear'',''extrap'')',SOCstr,OCVstr);
  cellParams.function.const.Uocv = str2func(OCVfn);
end

% Make electrode-level OCP function from MSMR parameter values
function Uocp = makeOCP(reg)
  F = 96485.3365; % Faraday's constant
  R = 8.3144621;  % Universal gas constant
  SOCvec = 0.001:0.001:0.999; % span entire electrode range

  U0       = reg.U0();
  X        = reg.X();
  omega    = reg.omega();
  
  % First, find relationship at Tref = 298.15K = 25 degC
  T = 298.15; f = F/(R*T);
  cost = @(U,Qdes) Qdes - sum(X./(1+exp(f*(U-U0)./omega)));
  OCVvec1 = 0*SOCvec;
  x0 = [0 5];
  for k = 1:length(OCVvec1)
    try
      OCVvec1(k) = fzero(cost,x0,[],SOCvec(k)); x0 = OCVvec1(k);
    catch
      OCVvec1(k) = NaN;
    end
  end
  
  % Next, find relationship at Tref = 299.15K = 26 degC
  T = 299.15; f = F/(R*T);
  cost = @(U,Qdes) Qdes - sum(X./(1+exp(f*(U-U0)./omega)));
  OCVvec2 = 0*SOCvec;
  x0 = [0 5];
  for k = 1:length(OCVvec2)
    try
      OCVvec2(k) = fzero(cost,x0,[],SOCvec(k)); x0 = OCVvec2(k);
    catch
      OCVvec2(k) = NaN;
    end
  end
  
  SOCstr = sprintf(' %g;',SOCvec);   SOCstr = ['[',SOCstr(1:end-1),']'];
  OCVstr1 = sprintf(' %g;',OCVvec1); OCVstr1 = ['[',OCVstr1(1:end-1),']'];
  OCVstr2 = sprintf(' %g;',OCVvec2-OCVvec1); 
  OCVstr2 = ['[',OCVstr2(1:end-1),']'];
  
  UocpFn = sprintf('@(x,T) interp1(%s,%s+(T(:)-298.15).*%s,x,''linear'',''extrap'')', ...
                   SOCstr,OCVstr1,OCVstr2);
  Uocp = str2func(UocpFn);         
end

% Make electrode-level dOCP function from MSMR parameter values
function dUocp = makedOCP(reg)
  F = 96485.3365; % Faraday's constant
  R = 8.3144621;  % Universal gas constant
  SOCvec = 0.001:0.001:0.999; % span entire electrode range
  dUdQ = 0*SOCvec;

  U0       = reg.U0();
  X        = reg.X();
  omega    = reg.omega();
  
  % Find dUdQ using Uocp at Tref
  T = 298.15; f = F/(R*T);
  Uocp = reg.Uocp(SOCvec,T); % get Uocp at Tref
  for k = 1:length(Uocp)
    x = X./(1+exp(f*(Uocp(k)-U0)./omega)); % get x_j for this U
    dQdU = -F/R*sum(x.*(1-x./X)./omega);
    dUdQ(k) = 1/dQdU;
  end
  
  SOCstr = sprintf(' %g;',SOCvec); SOCstr = ['[',SOCstr(1:end-1),']'];
  dOCVstr = sprintf(' %g;',dUdQ);  dOCVstr = ['[',dOCVstr(1:end-1),']'];
  
  % Make function, adding in temperature-dependence
  dUocpFn = sprintf('@(x,T) interp1(%s,T(:).*%s,x,''linear'',''extrap'')', ...
                    SOCstr,dOCVstr);
  dUocp = str2func(dUocpFn);   
end
