function [xmin,fmin,counteval,out] = gaussAdapt(fitfun,xstart,inopts)

global amp;

% Hacking for non-linear SDP solver
global C;
global counteval;

%-------------------------------------------------------------------------
% Implementation of the Gaussian Adaptation algorithm for design centering
% Black-box optimization and adaptive MCMC sampling
%
% Input:
% fitfun: Name of the fitness/target function as string or function handle
% xstart: initial candidate solution/sample
% inopts: option structure that determines internal strategy parameters (see code for details)
%
% Output:
% xmin: minimum candidate solution found by GaA (when using GaA as optimizer)
% fmin: fitness value of the xmin (when using GaA as optimizer)
% counteval: Number of function evaluations
% out: Output structure storing all relevant information (see code for
% details)
%
% When using this code please cite:
%
% C. L. Mueller and I. F. Sbalzarini. Gaussian Adaptation revisited - an
% entropic view on Covariance Matrix Adaptation.
% In Proc. EvoStar, volume 6024 of Lecture Notes in Computer Science,
% pages 432?441, Istanbul, Turkey, April 2010. Springer.
%
% C. L. Mueller and I. F. Sbalzarini. Gaussian Adaptation as a unifying
% framework for continuous black-box optimization and adaptive Monte Carlo sampling.
% In Proc. IEEE Congress on Evolutionary Computation (CEC), Barcelona,
% Spain, July 2010.
%
% Christian L. Mueller
% MOSAIC group, ETH Zurich, Switzerland
%
% Version 01/11
%
%-------------------------------------------------------------------------

% dimension of the problem
N = length(xstart);

% Options defaults: Stopping criteria % (value of stop flag)
defopts.StopFitness   = -Inf;    % stop if f(xmin) < stopfitness, minimization';
defopts.MaxIter       = 1e4*(N); % maximal number of iterations/function evaluations';
defopts.TolX          = 1e-12;   % restart if history of xvals smaller TolX';
defopts.TolFun        = 1e-9;    % restart if history of funvals smaller TolFun';
defopts.TolR          = 1e-9;    % restart if step size smaller TolR';
defopts.TolCon        = 1e-9;    % restart if threshold and fitness converge';
defopts.BoundActive   = 0;       % Flag for existence of bounds
defopts.BoundPenalty  = 0;       % Flag for use of penalty term outside bounds;
defopts.LBounds       = -Inf;    % lower bounds, scalar or Nx1-vector';
defopts.UBounds       = Inf;     % upper bounds, scalar or Nx1-vector';
defopts.A             = [];      % Constraint matrix Ax<=b;
defopts.E             = [];      % Constraint ellipsoid (x-cE)'*E*(x-cE) <=1;
defopts.cE            = [];      % Constraint ellispoid center (x-cE)'E(x-cE) <=1;
defopts.bVec          = [];      % constraint vector Ax<=b;
defopts.bRestart      = 1;       % Flag for restart activation;
defopts.ThreshRank    = 0;       % Flag for threshold based on ranks (Experimental)
defopts.PopMode       = 0;       % Flag for population mode (ToDo);
defopts.Display       = 'off';   % Display 2D landscape while running';
defopts.Plotting      = 'on';    % Plot progress while running';
defopts.VerboseModulo = 1e3;     % >=0, command line messages after every i-th iteration';
defopts.SavingModulo  = 1e2;     % >=0, saving after every i-th iteration';
defopts.bSaving       = 'on';    % [on|off] save data to file';
defopts.bSaveCov      = 0;       % save covariance matrices' ('1' results in huge files);
defopts.funArgs       = [];      % give additional target function arguments
%(scalar, vectors or matrix) if necessary

% Default options for algorithmic parameters
defopts.valP    = 1/exp(1);         % 0.37... Hitting probability
defopts.r       = 1;                % Step size of the initial covariance
defopts.initQ   = eye(N,N);         % Initial Cholesky matrix
defopts.MaxR    = Inf;              % maximal allowed step size
defopts.MinR    = 0;                % minimal allowed step size
defopts.MaxCond = 1e20*N;            % maximal allowed condition
defopts.N_mu    = exp(1)*N;         % Mean adaptation weight
defopts.N_C     = (N+1)^2/log(N+1); % Matrix adaptation weight
defopts.N_T     = exp(1)*N;         % Constraint adaptation weight
defopts.inc_T   = 2;                % Factor for N_T increase at restart (optimization)
defopts.beta    = 1/defopts.N_C;    % Step size increase/decrease factor
defopts.ss      = 1 + defopts.beta*(1-defopts.valP); % Expansion upon success
defopts.sf      = 1 - defopts.beta*(defopts.valP);   % Contraction otherwise

% Option for optimization/design-centering/MCMC sampling mode
% mode = 0 design centering
% mode = 1 optimization
% mode = 2 MCMC sampling

defopts.mode = 1;

% Option for initial threshold
% In the design centering mode c_T is constant. Values smaller than this value are
% considered as feasible points
% In the optimization mode this value is adapted
% In the MCMC mode this threshold is neglected
defopts.c_T = Inf;

% Option for initial sample point matrix
defopts.initX = [];

% ---------------------- Handling Input Parameters ----------------------

if isempty(fitfun)
    error('Objective function not determined');
end
if ~ischar(fitfun)
    error('first argument FUN must be a string');
end

if nargin < 2
    xstart = [];
end

if isempty(xstart)
    error('Initial search point, and problem dimension, not determined');
end

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    inopts = [];
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end


% ---------------------- Setup algorithmic Parameters ----------------------

% Options for algorithmic parameters
P = opts.valP;

if ~isfield(inopts,'N_T') && isfield(inopts,'N_C')
    opts.N_T = opts.N_C/2;
end

if ~isfield(inopts,'beta') && isfield(inopts,'N_C')
    opts.beta = 1/(opts.N_C);
end

if ~isfield(inopts,'ss')
    opts.ss = 1 + opts.beta*(1-opts.valP); % Expansion
end

if ~isfield(inopts,'sf')
    opts.sf = 1 - opts.beta*(opts.valP); % Contraction
end

if isfield(inopts,'UBounds') && isfield(inopts,'LBounds')
    opts.BoundActive  = 1;
    
    if ~isfield(inopts,'MaxR')
        opts.MaxR = 0.5*max(opts.UBounds-opts.LBounds);
    end
    
    % Make vector out of scalar input
    if size(opts.UBounds,1)==1
        opts.UBounds = ones(N,1)*opts.UBounds;
    end
    if size(opts.LBounds,1)==1
        opts.LBounds = ones(N,1)*opts.LBounds;
    end
    
end

% Algorithmic parameters

N_mu = opts.N_mu;
N_C = opts.N_C;
N_T = opts.N_T;

beta = opts.beta;
ss = opts.ss;
sf = opts.sf;

r = opts.r;
r_init = r;
rMax = opts.MaxR;
rMin = opts.MinR;

inc_T = opts.inc_T;

condMax = opts.MaxCond;

bAdapt=0;

P_emp = 0;

% Display the current options
% opts
% pause

% Initialize dynamic (internal) strategy parameters
Q = opts.initQ;
C = r^2*Q'*Q;
[eigVecs,eigVals] = eig(Q*Q');

C_old = C;


% Condition of initial Q
condC = cond(Q*Q');

% Initialize plotting figure
if strcmp(opts.Plotting,'on')
    
    figure(42)
    hold on
    grid on
    lastInd = 1;
    lastSave = 1;
    lastEigVals = ones(N,1);
    
end

% Compute function for 2D plotting
if strcmp(opts.Display,'on')
    
    if strcmp(fitfun,'benchmark_func')
        global initial_flag;
        initial_flag=0;
    end
    
    figure(23)
    
    try
        x = linspace(opts.LBounds(1),opts.UBounds(1),100);
        y = linspace(opts.LBounds(2),opts.UBounds(2),100);
        f = zeros(length(x),length(y));
    catch
        error('Displaying a function requires boundary setting')
    end
    
    i=1;
    j=1;
    for j=1:length(x)
        for i=1:length(y)
            if isempty(opts.funArgs)
                f(i,j)=feval(fitfun, [x(j);y(i)]);
            else
                f(i,j)=feval(fitfun, [x(j);y(i)],opts.funArgs);
            end
        end
    end
    
    [X,Y]=meshgrid(x,y);
    % meshc(Y,X,f)
    % hold on
    
    if opts.mode ~= 0
        
        if strcmp(fitfun,'fnoisysphere')
            contour(X,Y,f,10)
        else
            contour(X,Y,f,100)
        end
        
        if strcmp(fitfun,'benchmark_func')
            initial_flag=0;
        end
        
    else
        % Special display for design centering
        surf(X,Y,f)
        colormap('autumn')
        shading flat
        view(0,90)
        xlim([opts.LBounds(1) opts.UBounds(1)])
        ylim([opts.LBounds(2) opts.UBounds(2)])
    end
    
    hold on
    
end

% ----------------- Setup initial settings ---------------------------

% Determine start fitness value and threshold
if isempty(opts.funArgs)
    c_T=feval(fitfun, xstart);
else
    c_T=feval(fitfun, xstart,opts.funArgs);
end

arfitness=c_T;

% Design centering or optimization
if (opts.mode==0)
    if isfield(inopts,'c_T')
        c_T = opts.c_T;
        % Samples below threshold are considered feasible points
        threshold = c_T;
    end
end

% For testing issues and (possible) restarts save the initial values
arfitness_init = arfitness;
c_T_init = c_T;
N_mu_init = N_mu;
N_C_init = N_C;
N_T_init = N_T;

lastxAcc = xstart;
lastfAcc = arfitness_init;

lBnd=0;

% Linear constraints active
if ~(isempty(opts.A) && isempty(opts.bVec))
    lBnd=1;
    A = inopts.A;
    b = inopts.bVec;
    copts.thinInterval=1;
    copts.burnin=2*N;
    copts.bndEps = 1e-15;
end

eBnd=0;

if ~(isempty(opts.E) && isempty(opts.cE))
    eBnd=1;
    E = inopts.E;
    cE = inopts.cE;
    
    copts.thinInterval=1;
    copts.burnin=2*N;
    copts.bndEps = 1e-15;
end

% Transform lower/uppper bounds into matrix inequalities
bBnd = 0;
if opts.BoundActive && ~opts.BoundPenalty
    bBnd=1;
    
    if isfield(inopts,'UBounds') && isfield(inopts,'LBounds')
        A = [-eye(N,N);eye(N,N)];
        b = [-opts.LBounds;opts.UBounds];
    end
    
    if isfield(inopts,'UBounds') && ~isfield(inopts,'LBounds')
        A = eye(N,N);
        b = opts.UBounds;
    end
    
    if ~isfield(inopts,'UBounds') && isfield(inopts,'LBounds')
        A = -eye(N,N);
        b = -opts.LBounds;
    end
    
    copts.thinInterval=1;
    copts.burnin=2*N;
    copts.bndEps = 1e-12;
end

% Number of iterations
counteval = 1;

% Number of accepted points
numAcc = 0;

% Number of accepted points within a restart
rNumAcc = 0;

% ----------------- Setup output Parameters ---------------------------

mu = xstart;
mu_old = mu;
arx = xstart;

% Restart specific bestever
currxBest = xstart;
currfBest = arfitness;

% Overall bestever
xmin = xstart;
fmin = arfitness;
Qmin = Q;
rmin = r;

% History for rank-based selection (experimental)
histLen=10*N;%2*ceil(4+3*log(N));
histFAcc = zeros(histLen,1);
histXAcc = zeros(histLen,N);
histFAcc(:) = currfBest;
sortedHist = histFAcc;

if strcmp(opts.bSaving,'on')
    
    % Trace of raw samples and function values
    xRaw=zeros(opts.MaxIter/opts.SavingModulo,N);
    xRaw(1,:)=xstart;
    fRaw=zeros(opts.MaxIter/opts.SavingModulo,1);
    fRaw(1) = c_T;
    
    % Trace of accepted samples and function values
    xAcc=zeros(opts.MaxIter/opts.SavingModulo,N);
    xAcc(1,:)=xstart';
    fAcc=zeros(opts.MaxIter/opts.SavingModulo,1);
    fAcc(1) = arfitness;
    cntAcc=zeros(opts.MaxIter/opts.SavingModulo,1);
    
    % Trace of best samples and function values
    xBest=zeros(opts.MaxIter/opts.SavingModulo,N);
    xBest(1,:)=xstart;
    fBest=zeros(opts.MaxIter/opts.SavingModulo,1);
    fBest(1) = c_T;
    
    % Vector of step lengths
    rVec=zeros(opts.MaxIter/opts.SavingModulo,1);
    rVec(1)=r;
    
    % Vector of KL divergence
    KLVec=zeros(opts.MaxIter/opts.SavingModulo,1);
    KLVec(1)=0;
    KL=0;
    
    % Vector of mu
    muVec = zeros(opts.MaxIter/opts.SavingModulo,N);
    muVec(1,:)=mu;
    
    % Vector of acceptance thresholds
    c_TVec=zeros(opts.MaxIter/opts.SavingModulo,1);
    c_TVec(1)=c_T;
    
    % Vector of empirical acceptance probability
    P_empVec = zeros(opts.MaxIter/opts.SavingModulo,1);
    P_empVec(1)=opts.valP;
    
    % Cell array of Q matrices
    if opts.bSaveCov==1
        QCell = cell(opts.MaxIter/opts.SavingModulo,1);
        QCell{1} = Q;
    end
    
    % Vector of evaluation indices
    countVec = zeros(opts.MaxIter/opts.SavingModulo,1);
    countVec(1) = 1;
    
    % Cell array with stop flags
    stopFlagCell = cell(1,1);
    
end

% Final initialization
saveInd = 2;
stopFlag='';
cntInitX = 1;

% -------------------- Generation Loop --------------------------------

% Reserve the last evaluation for the fitness at the last mu
while counteval < (opts.MaxIter-1)
    
    arfitness_old=arfitness;
    C_old = C;
    counteval = counteval+1;
    
    if strcmp(opts.Display,'on') && (mod(counteval,opts.VerboseModulo)==0)
        figure(23)
        error_ellipse(C([1,N],[1,N]),mu([1,N],1),'style','k')
        grid on
    end
    
    
    % Handle boundaries
    penalty=0;
    
    if lBnd || bBnd
        arx_old = arx;
        %arx = mvntruncHR(mu,r,C./r^2,A,b,1);
        %arx = mvntruncLMI(mu,r,C./r^2,A,1);
        
        %arx = mvntrunc(mu,r,C./r^2,A,b,1,copts);
        arx = mvntruncFortran(mu,r,C./r^2,A,b,1,copts.burnin,copts.thinInterval);
        %arx = mvntruncFortran(mu,r,eigVecs,diag(eigVals),A,b,1,copts.burnin,copts.thinInterval);
        arz = Q'\((arx-mu)./r);
        
    elseif eBnd
        arx_old = arx;
        arx = mvntruncElli(mu,r,C,cE,E,1,copts);
        arz = Q'\((arx-mu)./r);
        
        %     elseif bBnd
        %
        %         arx_old = arx;
        %         arx = mvntrunc2(mu,C,[opts.LBounds,opts.UBounds],1,copts);
        %         arz = Q'\((arx-mu)./r);
        
    else
        
        % Generate a N(0,1) variable
        arz = randn(N,1);  % array of standard normally distributed mutation vectors
        
        % Sample from N(mu,C) distribution
        arx_old = arx;
        arx = mu + r * (Q * arz);
        
        % Box constraints active
        if opts.BoundActive
            arxNoBound=arx;
            arx=max([arx,opts.LBounds],[],2);
            arx=min([arx,opts.UBounds],[],2);
            if opts.BoundPenalty
                % Simple boundary handling
                penFac=mean(abs(sortedHist(1:N)));
                penalty=penFac*norm(arx-arxNoBound).^2;
            end
        end
        
    end
    % Plug-in local optimizer (experimental)
    % arx=fminunc(fitfun,arx);
    
    % Objective/target function call
    if isempty(opts.funArgs)
        arfitness=feval(fitfun, arx)+penalty;
        % Just for recording
        % muFitness=feval(fitfun, mu);
    else
        arfitness=feval(fitfun, arx,opts.funArgs)+penalty;
        % Just for recording
        % muFitness=feval(fitfun, mu,opts.funArgs);
    end
    
    % Save best sample and fitness values per restart
    if arfitness<currfBest
        currfBest = arfitness;
        currxBest = arx';
    end
    
    % Save best ever sample and fitness values across restarts
    if arfitness<fmin
        fmin = arfitness;
        xmin = arx;
        Qmin = Q;
        rmin = r;
    end
    
    if strcmp(opts.Plotting,'on') && ...
            (mod(counteval,opts.VerboseModulo)==0) && strcmp(opts.bSaving,'on')
        %tic
        % shift values such that the minimal value is 0
        if (opts.StopFitness  > -Inf)
            biasVal = -opts.StopFitness;
        else
            biasVal = 0;
        end
        
        currInd = counteval-1;
        currIndices = lastInd:opts.SavingModulo:currInd;
        currSaves = lastSave:saveInd-1;
        
        %[eigVecs,eigVals] = eig(C./r^2); % eigen decomposition for plotting
        eigValsD=diag(eigVals);
        
        try
            figure(42)
            subplot(2,2,1)
            plot([currIndices],muVec(currSaves,:)')
            title(['Current best: ',num2str(currfBest)])
            grid on
            hold on
            
            subplot(2,2,2)
            semilogy([currIndices],rVec(currSaves,:)')
            title(['Current step size: ',num2str(r),' and KL divergence: ',num2str(KL)])
            grid on
            hold on
            semilogy([currIndices],KLVec(currSaves,:)','r')
            
            
            subplot(2,2,3)
            semilogy([currIndices],biasVal+fBest([currSaves]))
            hold on
            grid on
            semilogy([currIndices],biasVal+c_TVec(currSaves,:)','r')
            title(['Current threshold: ',num2str(c_T)])
            
            % subplot(2,2,4)
            % plot([currIndices],P_empVec([currSaves]))
            % grid on
            % hold on
            % title(['Current acceptance prob: ',num2str(P_emp)])
            
            subplot(2,2,4)
            if N==2
                semilogy([lastInd,currInd],[1./sort(lastEigVals),1./sort(eigValsD)]','-')
            else
                semilogy([lastInd,currInd],[1./sort(lastEigVals),1./sort(eigValsD)],'-')
            end
            
            grid on
            hold on
            title(['Condition of C: ',num2str(max(eigValsD)/min(eigValsD))])
        end
        
        lastEigVals = eigValsD;
        
        lastInd = currInd;
        lastSave = saveInd-1;
        
        drawnow;
        %toc
    end
    
    % Design centering
    if (opts.mode==0)
        
        threshold = c_T;
        
        % Check whether to adapt or not
        bAdapt = (arfitness<threshold);
        
        % Optimization
    elseif (opts.mode==1)
        
        % Check whether fitness/rank-based selection is used
        if opts.ThreshRank % (experimental)
            threshold = sortedHist(floor(opts.valP*histLen));
        else
            threshold = c_T;
        end
        
        % Check whether to adapt or not
        bAdapt = (arfitness<threshold);
        
        % MCMC Sampling
    elseif (opts.mode==2)
        
        % Metropolis criterion for adaptive MCMC (M-GaA)
        accProb = min(1,arfitness/arfitness_old);
        % accProb = min(1,exp(arfitness_old-arfitness));
        bAdapt = (rand<accProb);
    end
    
    
    % Adapt moments upon acceptance (for all modes of operation)
    if bAdapt
        
        % Count accepted solutions
        numAcc = numAcc + 1;
        rNumAcc = rNumAcc + 1;
        
        % Expand r by factor ss
        r = r*ss;
        
        % Upper bound for step size to avoid numerical errors
        r = max(min(r,rMax),rMin);
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Adapt threshold c_T when optimizing
        %%%%%%%%%%%%%%%%%%%%%%
        if (opts.mode==1)
            
            % Update history of accepted f and x
            histFAcc(2:end)=histFAcc(1:end-1);
            histFAcc(1)=arfitness;
            sortedHist = sort(histFAcc,'ascend');
            
            histXAcc(2:end,:)=histXAcc(1:end-1,:);
            histXAcc(1,:)=arx';
            
            if opts.ThreshRank %(experimental)
                
                % threshold for last accepted sample
                c_T = threshold;
                
            else
                
                % Standard fitness-dependent threshold decrease
                c_T = (1-1/N_T)*c_T + arfitness/N_T;
                
                % Linear threshold decrease (experimental)
                % c_T = (1-counteval/opts.MaxIter)*300;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Adapt mu
        %%%%%%%%%%%%%%%%%%%%%%
        gamma = 1;%counteval;
        
        mu_old=mu;
        mu = (1-1/(gamma*N_mu))*mu + arx/(gamma*N_mu);
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Adapt covariance
        %%%%%%%%%%%%%%%%%%%%%%
        
        % This would be a direct update of C
        % C_direct = (1-1/N_C)*C + ((mu-arx)*(mu-arx)')/N_C
        
        
        % In order to avoid numerical errors
        if (condC<condMax) && (condC > 0)
            
            deltaC =  (1-1/N_C)*eye(N,N) + arz*arz'./N_C;
            % set deltaC to identity if covariance adaptation is not
            % required (experimental)
            % deltaC =  eye(N,N);
            
            % Adapt Q
            
            deltaC = triu(deltaC)+triu(deltaC,1)';  % enforce symmetry
            [B,eigValsD] = eig(deltaC);               % eigen decomposition, B==normalized eigenvectors
            deltaQ = B*diag(sqrt(diag(eigValsD)))*B'; % deltaQ contains standard deviations now
            
            if any(imag(deltaQ(:)))
                deltaQ
                error('deltaQ is imaginary')
            end
            
            % Update Cholesky matrix
            Q = Q*deltaQ;
            
            if any(imag(Q(:)))
                Q
                error('Q is imaginary')
            end
            
        else
            
            if (mod(counteval,opts.VerboseModulo)==0)
                disp( '-------------------------------------------');
                disp(' Condition of C is too large and regularized');
                disp( '-------------------------------------------');
            end
            
            % Regularize Q
            Q = Q + 1/N*eye(N,N);
            
        end
        
        % Normalize volume of Q to 1
        detQ=det(Q);
        Q = Q./((detQ).^(1/N));
        
        % Update condition of C
        [eigVecs,eigVals] = eig(Q*Q');
        condC = max(diag(eigVals))/min(diag(eigVals));
        
        if any(diag(eigVals)<0)
            saveInd
            xRaw(max(1,saveInd-10):saveInd,:)
            eigVecs
            eigVals
            condC
            detQ
            Q*Q'
            warning('C contains negative eigenvalues')
            eigVals(eigVals<0)=1e-3;
            Q = eigVecs'*eigVals*eigVecs;
        end
        
        
        % New full covariance matrix
        C = r^2*Q*Q';
        
        % GaA with vanishing adaptation (experimental)
        % N_C = counteval * N_C_init;
        
        lastxAcc = arx;
        lastfAcc = arfitness;
        
        
    else
        
        % Contract r by factor sf
        r=r*sf;
        
        % New full covariance matrix
        C = sf^2*C;
        
        
        % When MCMC sampling is on
        % set new position and fitness to old one (Metropolis criterion)
        if (opts.mode == 2)
            arfitness = arfitness_old;
            arx = arx_old;
        end
        
    end
    
    
    % Save accepted points (either really accepted or old one) TODO
    if strcmp(opts.bSaving,'on') && (mod(counteval,opts.SavingModulo)==0)
        cntAcc(saveInd,1)=counteval;
        xAcc(saveInd,:)=lastxAcc';
        fAcc(saveInd,1)=lastfAcc;
    end
    
    
    % Restart and stop flag checks in optimization mode
    if (opts.mode == 1)
        
        % Break, if fitness is good enough
        if arfitness <= opts.StopFitness
            fmin=arfitness;
            xmin=arx;
            stopFlag='fitness';
            break;
        end
        
        
        % Restart when history of accepted samples is converging
        
        %%%%% TolFun %%%%%
        currTolFun = abs(max(histFAcc)-min(histFAcc));
        if (currTolFun<opts.TolFun) && rNumAcc>histLen
            stopFlag='TolFun';
        end
        
        
        %%%%% TolX %%%%%
        currTolX = norm(histXAcc(1,:)-histXAcc(end,:));
        if (currTolX<opts.TolX)&& rNumAcc>histLen
            stopFlag='TolX';
        end
        
        %%%%% TolR %%%%%
        if (r<opts.TolR)
            stopFlag='TolR';
        end
        
        %%%%% TolCon %%%%%
        currTolCon = abs(c_T-currfBest);
        if N_T~=1 && currTolCon<opts.TolCon && rNumAcc>histLen
            stopFlag='TolCon';
        end
        
        % Restart when Restart flag is active AND stop flags are active
        if opts.bRestart && ~(strcmp(stopFlag,''))
            
            % Show reason for restart
            disp([num2str(counteval),': Restart due to: ',stopFlag]);
            
            % Record stop flags
            if size(stopFlagCell,1)==1
                stopFlagCell{1,1}=stopFlag;
            else
                stopFlagCell{end+1,1}=stopFlag;
            end
            
            % Reset stopFlag
            stopFlag='';
            
            % Reset Gaussian parameters
            r = r_init;
            Q = eye(N,N);
            C = r^2*Q'*Q;
            C_old = C;
            
            % Create new start point and handle boundaries
            if opts.BoundActive==1
                if ~isempty(opts.initX);
                    
                    cntInitX = max(1,mod(cntInitX,size(opts.initX,2)+1));
                    mu = opts.initX(:,cntInitX);
                    mu_old = mu;
                    arx = mu;
                    cntInitX = cntInitX+1;
                    disp('Use user-defined starting points')
                    
                elseif isfield(inopts,'UBounds') && isfield(inopts,'LBounds')
                    mu = opts.LBounds + rand(N,1).*(opts.UBounds-opts.LBounds);
                    arx = mu;
                    mu_old = mu;
                else
                    mu = arx;
                    arx = mu;
                    mu_old = mu;
                end
                
            else
                mu = xstart;
                arx = xstart;
                mu_old = mu;
            end
            
            % New start value
            if isempty(opts.funArgs)
                c_T=feval(fitfun, arx);
            else
                c_T=feval(fitfun, arx, opts.funArgs);
            end
            
            counteval = counteval+1;
            arfitness = c_T;
            rNumAcc = 0;
            
            % Restart with larger N_T (see CEC 2010 paper)
            N_T = inc_T*N_T;
            
            % Reset history for rank-based selection/tolerance checks
            histFAcc = arfitness*ones(histLen,1);
            sortedHist = histFAcc;
            
            currfBest=arfitness;
            histXAcc = repmat(arx',histLen,1);
            
        end
    end
    
    % Empirical acceptance probability
    P_emp = numAcc/counteval;
    
    % Save data
    if strcmp(opts.bSaving,'on') && (mod(counteval,opts.SavingModulo)==0)
        
        rVec(saveInd)=r;
        muVec(saveInd,:)=mu;
        c_TVec(saveInd)=c_T;
        P_empVec(saveInd)=P_emp;
        invC = inv(C);
        KL = 1/2*(log(det(C))-log(det(C_old)) + trace(invC*C_old) + (mu-mu_old)'*invC*(mu-mu_old) - N);
        KLVec(saveInd) = KL;
        
        % Cell array of Q matrices if covariances are stored
        if opts.bSaveCov==1
            QCell{saveInd}=Q;
        end
        
        fRaw(saveInd)=arfitness;
        xRaw(saveInd,:)=arx;
        
        fBest(saveInd)=currfBest;
        xBest(saveInd,:)=currxBest;
        
        countVec(saveInd)=counteval;
        
        saveInd = saveInd+1;
        
    end
    
    if (mod(counteval,opts.VerboseModulo)==0)
        disp( '-------------------------------------------');
        disp([' Number of iterations: ',num2str(counteval)]);
        disp([' P_acc:                ',num2str(P_emp)]);
        disp([' Search radius:        ',num2str(r)]);
        if opts.mode == 1
            disp([' Current fitness:      ',num2str(arfitness)]);
            disp([' Threshold:            ',num2str(c_T)]);
            disp([' TolCon:               ',num2str(currTolCon)]);
            disp([' TolX:                 ',num2str(currTolX)]);
        end
    end
    
    % No restart but stop flag is set
    if ~(strcmp(stopFlag,''))
        break;
    end
    
end % while, end generation loop

% Final sample is reserved for the final mu vector
if (counteval==(opts.MaxIter-1))
    
    % Final objective function call is reserved for mu
    if isempty(opts.funArgs)
        muFitness=feval(fitfun, mu);
    else
        muFitness=feval(fitfun, mu,opts.funArgs);
    end
    
    if muFitness <= opts.StopFitness
        currfBest = muFitness;
        currxBest = mu;
        stopFlag='fitness';
        if muFitness <=fmin
            fmin = muFitness;
            xmin = mu;
            Qmin = Q;
            rmin=r;
        end
    else
        stopFlag='MaxIter';
    end
    counteval = counteval+1;
end

if strcmp(opts.bSaving,'on')
    
    % Save the last evaluation
    rVec(saveInd) = r;
    muVec(saveInd,:) = mu;
    xRaw(saveInd,:) = arx;
    fRaw(saveInd,:) = arfitness;
    c_TVec(saveInd) = c_T;
    P_empVec(saveInd) = P_emp;
    KLVec(saveInd) = KL;
    
    if opts.bSaveCov==1
        QCell{saveInd}=Q;
    end
    
    fBest(saveInd)=currfBest;
    xBest(saveInd,:)=currxBest;
    countVec(saveInd)=counteval;
    
    out.xRaw = xRaw(1:saveInd,:);
    out.fRaw = fRaw(1:saveInd,:);
    out.cntAcc=cntAcc(1:saveInd-1,:);
    out.fAcc=fAcc(1:saveInd-1,:);
    out.xAcc=xAcc(1:saveInd-1,:);
    out.xBest = xBest(1:saveInd,:);
    out.fBest = fBest(1:saveInd,:);
    out.P_empVec = P_empVec(1:saveInd,:);
    out.rVec = rVec(1:saveInd,:);
    out.muVec = muVec(1:saveInd,:);
    out.c_TVec = c_TVec(1:saveInd,:);
    out.KLVec = KLVec(1:saveInd,:);
    
    out.bestever.x = xmin;
    out.bestever.f = fmin;
    out.bestever.Q = Qmin;
    out.bestever.r = rmin;
    
    if opts.bSaveCov==1
        out.QCell = QCell(1:saveInd,:);
    end
    out.countVec = countVec(1:saveInd);
    
    % Record all stop flags
    if isempty(stopFlagCell{1,1})
        stopFlagCell{1,1}=stopFlag;
    else
        stopFlagCell{end+1,1}=stopFlag;
    end
    out.stopFlagCell = stopFlagCell;
    
end

% General output
out.stopFlag = stopFlag;
out.P_emp = numAcc/counteval;
out.lastQ = Q;
out.lastR = r;
out.lastMu = mu;
out.opts = opts;

% -------------------- Ending Message ----------------------------

disp(['          Acceptance probability: ' num2str(numAcc/counteval)]);
if (opts.mode == 1)
    disp([num2str(counteval) ' Current fitness: ' num2str(arfitness)]);
    disp(['          Best ever fitness: ' num2str(fmin)]);
end
disp(['          Stop flag: ' stopFlag]);

if strcmp(opts.Display,'on')
    figure(23)
    error_ellipse('C',C(1:2,1:2),'mu',mu(1:2),'style','b')
    hold on
    plot(mu(1,1),mu(2,1),'b.', 'MarkerSize',10)
end

% ------------------- Optimization test functions ----------------
% PARTLY TAKEN FROM Niko Hansen's CMA-ES code
% ----------------------------------------------------------------
function f=flin(x,c)
%global mu
%global A
try
    %H = dikinElli(A,mu,1);
    %xT=mu+chol(H)'*(x-mu);
    f=c'*x;
catch
    f=c'*x;
end


function f=fsphere(x)
f=sum(x.^2);

function f=flogsphere(x)
f=log(1+sum(x.^2));

function f=fschwefel(x)
f = sum(cumsum(x).^2);

function f=fcigar(x)
f = x(1)^2 + 1e6*sum(x(2:end).^2);

function f=fcigtab(x)
f = x(1)^2 + 1e8*x(end)^2 + 1e4*sum(x(2:(end-1)).^2);

function f=ftablet(x)
f = 1e6*x(1)^2 + sum(x(2:end).^2);

function f=felli(x)
N = size(x,1); if N < 2 error('dimension must be greater one'); end
f=1e6.^((0:N-1)/(N-1)) * x.^2;

function f=felli100(x)
N = size(x,1); if N < 2 error('dimension must be greater one'); end
f=1e4.^((0:N-1)/(N-1)) * x.^2;

function f=fplane(x)
f=x(1);

function f=frosen(x)
if size(x,1) < 2 error('dimension must be greater one'); end
f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);

function f=foverton(x)
if size(x,1) ~= 2 error('dimension must be two'); end
w=8;
f = w*abs(x(1,:).^2-x(2,:)) + (1-x(1,:)).^2;

function f=frastrigin(x)
N = size(x,1); if N < 2 error('dimension must be greater one'); end
amp=10;
f = amp * N + sum(x .^2 - amp * cos(2 * pi .* x),1);

function f=flograstrigin(x)
% Check for scale-invariance
N = size(x,1); if N < 2 error('dimension must be greater one'); end
%amp=10;
f = amp * N + sum(x .^2 - amp * cos(2 * pi .* x),1);
f=log(1+f);

function f=fcos(x,amp)
N = size(x,1); if N < 2 error('dimension must be greater one'); end
f = amp * N - amp * sum(cos(2 * pi .* x),1);

function f=fdsphere(x,params)
% Lunacek et al. 2008
s=params(1);
d=params(2);
N=length(x);
mu1=2.5*ones(N,1);
mu2=-sqrt((mu1.^2-d)./s);
f=min(sum((x-mu1).^2),d*N+s*sum((x-mu2).^2));

function f=fdrastrigin(x,params)
% Lunacek et al. 2008
N = length(x);
mu1 = 2.5*ones(N,1);

if length(params)==2
    s=params(1);
    d=params(2);
    mu2=-sqrt((mu1.^2-d)./s);
else
    s=params(1);
    mu2 = -2.5*ones(N,1);
    d = mu1(1).^2-s*mu2(1).^2;
end

f1=min(sum((x-mu1).^2),d*N+s*sum((x-mu2).^2));
x_mu=(x-mu1);
f2= 10*sum(1-cos(2 * pi .*x_mu),1);
f=f1+f2;

function f=fnoisysphere(x,sigma_eps)
% Additive noise
randNum = sigma_eps*randn;
f=norm(x)^2+randNum;
% Multiplicative noise
% randNum = -1.5+3*rand;
% f=norm(x)^2*(1+randNum);

function f=frand(x)
f=rand;

function f = MullerBrown(x)

[m,n] = size(x);
if m~=2
    error('Dimensionality must be 2')
end

MBMat = [-200,-1,0,-10,1,0;...
    -100,-1,0,-10,0,0.5;...
    -175 -6.5,11,-6.5,-0.5,1.5;...
    15,0.7,0.6,0.7,-1.0,1.0]';

f = zeros(n,1);

for i=1:n
    
    xVec = x(1,i).*ones(1,4);
    yVec = x(2,i).*ones(1,4);
    f(i) = sum(MBMat(1,:).*exp(MBMat(2,:).*(xVec-MBMat(5,:)).^2+MBMat(3,:).*...
        (xVec-MBMat(5,:)).*(yVec-MBMat(6,:))+MBMat(4,:).*(yVec-MBMat(6,:)).^2));
end



%%%%%%%%% Examples from Kjellstrom 1981-1999 %%%%%%%%%

function f=fkjellstrom1(x)
N = size(x,1); if N < 2 error('dimension must be greater one'); end
f = 20.0 - sum((x+0.5).^2) + exp(25*(sum(x.^2)-8));

function f=fkjellstrom2(x)
N = size(x,1); if N < 2 error('dimension must be greater one'); end
h = 0.01*(cos(x+1.982) + cos(2*x+5.720) + cos(3*x + 1.621) + cos(4*x + 0.823) + cos(5*x + 3.22));
f = prod(1+h(:));

% ------------------------------------------------------------------------
%%%%%%%%% Examples for test target distributions %%%%%%%%%
% ------------------------------------------------------------------------

function f=fHaario(x,type)
% from Haario et al. Comp Stat. 1999
global HaarioMat;
N = size(x,1);
C=eye(N,N);
switch type
    case 1
        C(1,1)=100;
    case 2
        C(1,1)=100;
        % TODO rotation
    case 3
        C(1,1)=100;
        b=0.03;
        x(2,1)=x(2,1)+b*x(1,1)^2-100*b;
    case 4
        C(1,1)=100;
        b=0.1;
        x(2,1)=x(2,1)+b*x(1,1)^2-100*b;
    otherwise
end
f = -mvnpdf(x,zeros(N,1),C);

function f=fLiangEx1(x)
% taken from Liang:2001 a 20-funnel distribution
global LiangExMat;
sigma=0.1;
w=0.05;
if length(x)~=2
    error('Wrong dimension')
else
    f=0;
    for i=1:2:39
        currCov = (x-LiangExMat(i:i+1))'*(x-LiangExMat(i:i+1));
        f=f+w*exp(-1/(2*sigma^2)*currCov);
    end
    f=-1/sqrt(2*pi*sigma)*f;
end

function f=fLiangEx2(x)
% taken from Liang:2001 a double-funnel distribution
N = length(x);
f = -1/3*mvnpdf(x,zeros(N,1),eye(N,N));
f = f - 2/3*mvnpdf(x,5*ones(N,1),eye(N,N));



function f=fNealsFunnel(x)
% from Neal's Slice sampler paper 2003
N = length(x);
%if N~=10
%    ('Input must be of dimension 10')
%end

v = x(1);
f = pdf('Normal',v,0,3);

for i=2:N
    f = f + pdf('Normal',x(i),0,sqrt(exp(v)));
end

function f=fSkilling(x)
% from Neal's Slice sampler paper 2003
N = length(x);
%if N~=20
%    ('Input must be of dimension 20')
%end

f1 = mvnpdf(x,zeros(N,1),0.1^2*eye(N,N));
f2 = mvnpdf(x,0.031*ones(N,1),0.01^2*eye(N,N));

f=f1+100*f2;



% ------------------------------------------------------------------------
%%% Example for design centering of linear/ellipsoidal feasible region %%%
% ------------------------------------------------------------------------


function f=fConRegion(x,params)
n = length(x);
A = params(:,1:n);
b = params(:,n+1);
if all(A*x+b<0) && all(x>=0)
    f=0;
else
    f=1;
end

function f=fConElli(x,params)
n = length(x);
A = params(1:n,1:n);
b = params(1:n,n+1);
if (x-b)'*A*(x-b)<1
    f=0;
else
    f=1;
end

function f=fConBox(x,params)
n = length(x);
sideLength = params(1);

if all(x>0 & x<sideLength)
    f=0;
else
    f=1;
end

function f=fConNSphere(x)
n = size(x,1);
A = eye(n,n);
b = zeros(n,3);
b(1,1)=-2;
b(1,2)=0;
b(1,3)=2;
if (x-b(:,1))'*A*(x-b(:,1))<1
    f=0;
elseif (x-b(:,2))'*A*(x-b(:,2))<1
    f=0;
elseif (x-b(:,3))'*A*(x-b(:,3))<1
    f=0;
else
    f=1;
end


% ---------------------------------------------------------------
% FUNCTIONS BELOW ARE TAKEN FROM Niko Hansen's CMA-ES code
% ---------------------------------------------------------------
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%


if nargin < 2 || isempty(defopts) % no default options available
    opts=inopts;
    return;
elseif isempty(inopts) % empty inopts invoke default options
    opts = defopts;
    return;
elseif ~isstruct(defopts) % handle a single option value
    if isempty(inopts)
        opts = defopts;
    elseif ~isstruct(inopts)
        opts = inopts;
    else
        error('Input options are a struct, while default options are not');
    end
    return;
elseif ~isstruct(inopts) % no valid input options
    error('The options need to be a struct or empty');
end

opts = defopts; % start from defopts
% if necessary overwrite opts fields by inopts values
defnames = fieldnames(defopts);
idxmatched = []; % indices of defopts that already matched
for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
    if isoctave
        for i = 1:size(defnames, 1)
            idx(i) = strncmp(lower(defnames(i)), lower(name), length(name));
        end
    else
        idx = strncmp(lower(defnames), lower(name), length(name));
    end
    if sum(idx) > 1
        error(['option "' name '" is not an unambigous abbreviation. ' ...
            'Use opts=RMFIELD(opts, ''' name, ...
            ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
        defname  = defnames{find(idx)};
        if ismember(find(idx), idxmatched)
            error(['input options match more than ones with "' ...
                defname '". ' ...
                'Use opts=RMFIELD(opts, ''' name, ...
                ''') to remove the field from the struct.']);
        end
        idxmatched = [idxmatched find(idx)];
        val = getfield(inopts, name);
        % next line can replace previous line from MATLAB version 6.5.0 on and in octave
        % val = inopts.(name);
        if isstruct(val) % valid syntax only from version 6.5.0
            opts = setfield(opts, defname, ...
                getoptions(val, getfield(defopts, defname)));
        elseif isstruct(getfield(defopts, defname))
            % next three lines can replace previous three lines from MATLAB
            % version 6.5.0 on
            %   opts.(defname) = ...
            %      getoptions(val, defopts.(defname));
            % elseif isstruct(defopts.(defname))
            warning(['option "' name '" disregarded (must be struct)']);
        elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
            opts = setfield(opts, defnames{find(idx)}, val);
            % next line can replace previous line from MATLAB version 6.5.0 on
            % opts.(defname) = inopts.(name);
        end
    else
        warning(['option "' name '" disregarded (unknown field name)']);
    end
end


% ---------------------------------------------------------------
% ---------------------------------------------------------------
function res = isoctave
% any hack to find out whether we are running octave
s = version;
res = 0;
if exist('fflush', 'builtin') && eval(s(1)) < 7
    res = 1;
end


