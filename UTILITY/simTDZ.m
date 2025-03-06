function simData = simTDZ(cellModel,freq,socPct,TdegC,varargin)
%SIMEIS Run time-domain EIS simulation in COMSOL over specified frequencies.
%
% data = SIMTDZ(cellModel,freq,socPct) simulates the full-order electrochemical
%   impedance spectroscopy (EIS) response a cell in COMSOL over the 
%   frequencies provided in the vector FREQ. CELLMODEL is a model of the 
%   cell obtained from LOADCELLPARAMS. SOCPCT is the state-of-charge of the
%   cell in percent [%]. DATA is a structure containing the time-domain
%   simulation data (see description below).
%
% data = SIMTDZ(cellModel,freq,socPct,TdegC) also specifies the temperature
%   at which to run the simulation. Default TdegC=25.
%
% data = SIMTDZ(...,'Vars',varsStruct) specifies the electrochemical variables 
%   to save in addition to Vcell and Iapp as well as the x-locations at which
%   to evalulate those variables. varsStruct is a structure containing
%   fields for each variable. The value of a field specifies the
%   x-locations at which to evalulate the specified variable.
%   Use the string value of 'mesh' to save all available x-locations.
%   Default varsStruct=struct('Phise',[0 3],'Phie',3).
%   NOTE: this function does not interpolate; if an x-location is specified
%   that is not in the COMSOL mesh, the closest x-location to that provided
%   is fetched (see which x-locations were fetched by accessing the `xlocs`
%   field of the returned DATA structure).
%
% data = SIMTDZ(...,'I',I) specifies the amplitude of the input sinusoid 
%   as I [A]. Default I=C/10.
%
% data = SIMTDZ(...,'Ns',Ns) specifies the sampling rate as the number of 
%   samples per period of the sinusoid [Sa/period]. Default Ns=32.
%
% data = SIMTDZ(...,'Nt',Nt) specifies the length of the transient interval 
%   in periods of the applied sinusoid [periods]. Default Nt=20.
%
% data = SIMTDZ(...,'Nss',Nss) specifies the length of the steady-state 
%   interval in periods of the applied sinusoid [periods]. Default Nss=8.
%
% data = SIMTDZ(...,'Verbose',true) outputs status messages to the command
%   window.
%
% The output, DATA, is a structure with the following fields:
%    ss = structure array containing steady-state waveforms. ss(k) is
%      a structure with the following fields corresponding to the kth
%      frequency in the input frequency vector freq:
%      - time: simulation time
%      - Vcell: vector of cell voltage vs time
%      - Iapp: vector of applied current vs time
%      - {Vars}: fields for each variable specified in the 'Vars' input
%        cell array. The value is a matrix whose first dimension (rows)
%        correspond to time and whose second dimension corresponds to
%        x-locations (columns).
%    xlocs = structure mapping variables to x-locations in the cell
%      sandwich.
%    param = structure of parameters supplied to the function.
%
% -- Changelog --
% 2025.02.21 | Updates for CDC paper | WH
% 2023.04.04 | Created | Wesley Hileman <whileman@uccs.edu>

% Fetch and validate input arguments.
isinteger = @(x)floor(x)==x;
parser = inputParser;
parser.addRequired('cellModel',@(x)isscalar(x));
parser.addRequired('freq',@(x)isvector(x)&&all(x>0));
parser.addRequired('socPct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addOptional('TdegC',25,@(x)isscalar(x));
parser.addParameter('Vars',struct('Phise',[0 3],'Phie',3),@(x)isstruct(x));
parser.addParameter('I',[],@(x)isscalar(x)&&x>0);
parser.addParameter('Ns',32,@(x)isscalar(x)&&isinteger(x)&&x>0);
parser.addParameter('Nt',20,@(x)isscalar(x)&&isinteger(x)&&x>0);
parser.addParameter('Nss',8,@(x)isscalar(x)&&isinteger(x)&&x>0);
parser.addParameter('Verbose',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('OptSimFOM',struct,@isstruct);
parser.parse(cellModel,freq,socPct,TdegC,varargin{:});
p = parser.Results;  % structure of validated arguments

if isempty(p.I)
    % Use C/10 rate as default.
    p.I = cellModel.function.const.Q()/10;
end

if p.Verbose
    fprintf('Starting time-domain EIS simulation\n');
end
FOM = genFOM(cellModel);

% Run simulation at each frequency (backwards to avoid need to pre-allocate
% structure arrays).
ff = p.freq(:);
varnames = fieldnames(p.Vars);
ssdata = [];
xlocsStruct = struct;
for k = length(ff):-1:1
    f0 = ff(k);

    % Generate time vector.
    kk = [0:p.Nt*p.Ns,(p.Nt*p.Ns+1):(p.Nt+p.Nss)*p.Ns]; % Discrete-time vector [sample number].      
    fs = p.Ns*f0; % Sampling rate [Sa/s].
    tt = kk/fs; % Time vector [s].
    ss = kk>p.Nt*p.Ns; % Logical indicies to steady-state interval.

    % Perform simulation in COMSOL.
    simspec.time = tt;
    simspec.mag = p.I;
    simspec.freq = f0;
    simspec.SOC0 = p.socPct;
    simspec.T = p.TdegC;
    simspec.TSHIFT = 0; % no need to shift for initial discontinuity (sine)
    if p.Verbose
        fprintf('Running f=%8.3gHz (%d of %d)...\n', ...
            f0,length(ff)-k+1,length(ff));
    end
    [~,sim] = simFOM(FOM,simspec,'InputType','sin',p.OptSimFOM);
    if p.Verbose
        fprintf('done!\n');
    end

    % Save steady-state waveforms.
    ssdata(k).param.f0 = f0;
    ssdata(k).param.fs = fs;
    ssdata(k).param.N = length(tt(ss));
    ssdata(k).time = tt(ss);
    ssdata(k).time = ssdata(k).time - ssdata(k).time(1); % start time at zero
    ssdata(k).Iapp = sim.Iapp(ss);
    ssdata(k).Vcell = sim.Vcell(ss);
    for i = 1:length(varnames)
        varname = varnames{i};
        timePositionData = sim.(varname);
        xlocs = sim.xLocs.(varname);      % COMSOL mesh locations
        xlocsDesired = p.Vars.(varname);  % desired x-locations
        if (ischar(xlocsDesired)||isstring(xlocsDesired))&&strcmpi(xlocsDesired,'mesh')
            % Store variable at all mesh locations.
            ssdata(k).(varname) = timePositionData(ss,:);
        else
            [~,indxx] = min(abs(xlocs(:)-xlocsDesired(:)'));
            ssdata(k).(varname) = timePositionData(ss,indxx);
            xlocs = xlocs(indxx); % true x-locations
        end
        if ~isfield(xlocsStruct,varname)
            xlocsStruct.(varname) = xlocs;
        end
    end
end

if p.Verbose
    fprintf('Finished time-domain EIS simulation\n');
end

% Collect results.
simData.ss = ssdata;
simData.xlocs = xlocsStruct;
simData.param = p;
simData.origin__ = 'simTDZ';

end