% function ROMout = simROMpack(ROMs,simData,architecture,method)
% 
% Inputs:
%   ROM          = matrix of reduced-order models created with an XRA
%   simData      = simulation profile loaded using "loadInput"
%   architecture = either 'SCM' or 'PCM'; the pack configuraion
%   method       = in the future, switch between output/model/non blending
%                  ("outBlend", "mdlBlend","nonBlend"). Presently, only
%                  "outBlend" is implemented.
%
% Output:
%   ROMout       = matrix of data structures with results from simulation
%
% This utility function executes a Ns by Np battery pack using
% reduced-order model and the user-prefered blending method
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

function ROMout = simROMpack(ROMs,simData,architecture,method)

  % Update progress
  fprintf('----------------------------------------------------------\n');
  fprintf('Start to simulate ROM-based pack using %s (%s)\n',architecture,datestr(now));

  switch upper(method)
    case 'OUTBLEND'      
      ROMout = outBlendPack(ROMs,simData,architecture);
    otherwise 
      error('simROMpack presently works only with output blending.')
  end

  fprintf('Finished: %s\n',datestr(now));
  fprintf('----------------------------------------------------------\n');
end