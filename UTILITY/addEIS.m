function model = addEIS(model)

  %% Set perturbation level... actual value should not matter because ...
  model.param.set('Iper', '1 [A]'); % COMSOL performs linear study
  
  %% Redefine temp
  model.param.set('T', '298.15 [K]', 'Ambient temperature');
  model.variable('varPos3').set('Vcell', 'phi_s-Iper*RcFN(1,T)');  

  %% Enable automatic equation forms to accommodate frequency inputs
  model.physics('phi_s').prop('EquationForm').setIndex('form', 'Automatic', 0);
  model.physics('phi_e').prop('EquationForm').setIndex('form', 'Automatic', 0);
  model.physics('theta_e').prop('EquationForm').setIndex('form', 'Automatic', 0);
  model.physics('if').prop('EquationForm').setIndex('form', 'Automatic', 0);
  model.physics('idl').prop('EquationForm').setIndex('form', 'Automatic', 0);
  model.component('mod2d').physics('theta_s').prop('EquationForm').setIndex('form', 'Automatic', 0);

  %% Add perturbation input to phi_s
  model.physics('phi_s').create('flux2', 'FluxBoundary', 0);
  model.physics('phi_s').feature('flux2').label('Flux/Source 2 - For perturbation study');
  model.physics('phi_s').feature('flux2').setIndex('g', 'linper(-Iper*xnorm)', 0);
  model.physics('phi_s').feature('flux2').selection.set(4);
  
  %% Create study 2, the perturbation study
  model.study.create('std2');
  model.study('std2').create('param', 'Parametric');
  model.study('std2').create('frlin', 'Frequencylinearized');
  model.study('std2').feature('frlin').set('useadvanceddisable', true);
  model.study('std2').feature('frlin').set('disabledphysics', {'phi_s/flux1'}); 
  model.study('std2').label('Study 2 - Frequency domain perturbation');
  model.study('std2').feature('param').set('pname', {'z0'});
  model.study('std2').feature('param').set('plistarr', {'0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1'});
  model.study('std2').feature('param').set('punit', {''});
  model.study('std2').feature('frlin').set('plist', '10^{range(log10(10000),-(1/10),log10(0.001))}');   
  
  %% Create solution 2 from the perturbation study    
  model.sol.create('sol2');
  model.sol('sol2').study('std2');
  model.sol('sol2').attach('std2');
  model.sol('sol2').create('st1', 'StudyStep');
  model.sol('sol2').create('v1', 'Variables');
  model.sol('sol2').create('s1', 'Stationary');
  model.sol('sol2').feature('s1').create('p1', 'Parametric');
  model.sol('sol2').feature('s1').create('fc1', 'FullyCoupled');
  model.sol('sol2').feature('s1').feature.remove('fcDef');
  model.sol('sol2').feature('st1').label('Compile Equations: Frequency-Domain Perturbation');
  model.sol('sol2').feature('v1').label('Dependent Variables 1.1');
  model.sol('sol2').feature('v1').set('clistctrl', {'p1'});
  model.sol('sol2').feature('v1').set('cname', {'freq'});
  model.sol('sol2').feature('v1').set('clist', {'10^{range(log10(10000),-(1/10),log10(0.001))}[Hz]'});
  model.sol('sol2').feature('s1').label('Stationary Solver 1.1');
  model.sol('sol2').feature('s1').set('nonlin', 'linper');
  model.sol('sol2').feature('s1').set('storelinpoint', true);
  model.sol('sol2').feature('s1').set('probesel', 'none');
  model.sol('sol2').feature('s1').feature('dDef').label('Direct 1');
  model.sol('sol2').feature('s1').feature('aDef').label('Advanced 1');
  model.sol('sol2').feature('s1').feature('p1').label('Parametric 1.1');
  model.sol('sol2').feature('s1').feature('p1').set('pname', {'freq'});
  model.sol('sol2').feature('s1').feature('p1').set('plistarr', {'10^{range(log10(10000),-(1/10),log10(0.001))}'});
  model.sol('sol2').feature('s1').feature('p1').set('punit', {'Hz'});
  model.sol('sol2').feature('s1').feature('p1').set('pcontinuationmode', 'no');
  model.sol('sol2').feature('s1').feature('p1').set('uselsqdata', false);
  model.sol('sol2').feature('s1').feature('fc1').label('Fully Coupled 1.1');
  model.sol('sol2').runAll;
  
  model.sol.create('sol4');
  model.sol('sol4').study('std2');
  model.sol('sol4').label('Parametric Solutions 1');
  
  model.batch.create('p1', 'Parametric');
  model.batch('p1').create('so1', 'Solutionseq');
  model.batch('p1').study('std2');
  model.batch('p1').set('control', 'param');
  model.batch('p1').set('pname', {'z0'});
  model.batch('p1').set('plistarr', {'0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1'});
  model.batch('p1').set('punit', {''});
  model.batch('p1').set('err', true);
  model.batch('p1').feature('so1').set('seq', 'sol2');
  model.batch('p1').feature('so1').set('psol', 'sol4');
  model.batch('p1').feature('so1').set('param', {'"z0","0.1"' '"z0","0.2"' '"z0","0.3"' '"z0","0.4"' '"z0","0.5"' '"z0","0.6"' '"z0","0.7"' '"z0","0.8"' '"z0","0.9"' '"z0","1"'});
  model.batch('p1').attach('std2');
  model.batch('p1').run;
  
  %% Make Nyquist plot of cell impedance
  model.result.create('pg39', 'PlotGroup1D');
  model.result('pg39').set('data', 'dset5');
  model.result('pg39').create('ptgr1', 'PointGraph');
  model.result('pg39').feature('ptgr1').set('data', 'dset5');
  model.result('pg39').feature('ptgr1').selection.set(4);
  model.result('pg39').feature('ptgr1').set('expr', '-imag(lindev(Vcell)/lindev(linper(-Iper)))');
  model.result('pg39').label('Impedance Nyquist Plot');
  model.result('pg39').set('looplevelinput', {'all' 'manual'});
  model.result('pg39').set('looplevel', {'1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71' '1'});
  model.result('pg39').set('xlabel', ['real(Z) (' 'ohm' ')']);
  model.result('pg39').set('ylabel', ['-imag(Z) (' 'ohm' ')']);
  model.result('pg39').set('xlabelactive', false);
  model.result('pg39').set('ylabelactive', false);
  model.result('pg39').feature('ptgr1').label('Impedance at different SOC');
  model.result('pg39').feature('ptgr1').set('xdata', 'expr');
  model.result('pg39').feature('ptgr1').set('xdataexpr', 'real(lindev(Vcell)/lindev(linper(-Iper)))+linpoint(RcFN(3,T))');
  model.result('pg39').feature('ptgr1').set('xdataunit','ohm');
  
  %% Make text export of Nyquist-plot data
  model.result.export.create('plot3', 'Plot');
  model.result.export('plot3').label('Export_Nyquist');
  model.result.export('plot3').set('plotgroup', 'pg39');
  model.result.export('plot3').set('plot', 'ptgr1');
  model.result.export('plot3').set('filename', 'impedance.txt');
end