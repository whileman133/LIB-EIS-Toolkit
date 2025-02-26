function [FOMout,model] = simEIS(model,SOCs,freqs)
  fprintf('Running study...\n');
  
  freqstr = sprintf('%g ',freqs);
  model.study('std2').feature('frlin').set('plist', freqstr);   
  model.study('std2').feature('param').set('plistarr', {num2str(SOCs)});
  model.batch('p1').set('plistarr', {num2str(SOCs)});

  z0str = {regexp(sprintf('"z0","%d"#',SOCs),'#','split')};
  model.batch('p1').feature('so1').set('param', z0str(1:end-1));
  model.batch('p1').attach('std2');
  model.batch('p1').run;

  model.sol('sol2').feature('v1').set('clist', freqstr);
  model.sol('sol2').feature('s1').feature('p1').set('plistarr', freqstr);
  
  model.sol('sol2').runAll;

  %% Gather output
  FOMout = [];
  FOMout.z0 = SOCs;
  FOMout.freqs = freqs;

  runData = [];
  for theZ = 1:length(SOCs)    
    out = mpheval(model,'lindev(Vcell)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    data_V = out.d1(:,end);
    out = mpheval(model,'lindev(linper(-Iper))','dataset','dset5','Edim',0,'outersolnum',theZ);      
    data_I = out.d1(:,end);
    out = mpheval(model,'RcFN(1,T)','dataset','dset5','Edim',0);
    Rc = out.d1(end);
    Z = data_V./data_I + Rc;

    phi_s    = mpheval(model,'lindev(phi_s)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    phi_e    = mpheval(model,'lindev(phi_e)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    phi_se   = mpheval(model,'lindev(phi_s-phi_e)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    theta_ss = mpheval(model,'lindev(thetass)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    theta_e  = mpheval(model,'lindev(theta_e)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    i_f      = mpheval(model,'lindev(if)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    i_dl     = mpheval(model,'lindev(idl)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    % Note: If I query "ifdl" instead of "if+idl", the results are poor
    % near separator boundaries(?)
    i_fdl    = mpheval(model,'lindev(if+idl)','dataset','dset5','Edim',0,'outersolnum',theZ);      
    eta      = mpheval(model,'lindev(eta)','dataset','dset5','Edim',0,'outersolnum',theZ);      

    runData(theZ).data_V = data_V; %#ok<*AGROW>
    runData(theZ).data_I = data_I;
    runData(theZ).Rc     = Rc;
    runData(theZ).Z      = Z;

    runData(theZ).phi_s    = phi_s;
    runData(theZ).phi_e    = phi_e;
    runData(theZ).phi_se   = phi_se;
    runData(theZ).theta_ss = theta_ss;
    runData(theZ).theta_e  = theta_e;
    runData(theZ).i_f      = i_f;
    runData(theZ).i_dl     = i_dl;
    runData(theZ).i_fdl    = i_fdl;
    runData(theZ).eta      = eta;    
  end
  FOMout.runData = runData;
end
