clear inopts
close all

% Function
funName = 'fmisteli'

N = 2;
numIter = 1;
MaxIter = N*1.5e2;

% number of iterations
for j=1:numIter
    
    disp(['Dimenison: ',num2str(N), ', Iteration: ',num2str(j)])
    
    % Set-up of GaA options
    inopts.VerboseModulo=1e1;
    inopts.bSaving='on';
    inopts.SavingModulo=1;
    
    inopts.r = 1;
    inopts.StopFitness = 1e-6;
    
    inopts.MaxIter=MaxIter;
    inopts.Display = 'off';
    inopts.Plotting = 'on';
    inopts.valP = 1/exp(1);
    
    %inopts.N_mu = 1;
    %inopts.N_T = 1;
    
    % Start value
    xstart = ones(N,1)
    [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
    
end
