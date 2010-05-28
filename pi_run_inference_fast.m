function [Phat,F,V] = pi_run_inference_fast(fast,F,V,options)

%% set up variables
V.est_n=1;          % estimate spiking parameters
V.est_h=1;          % estimate spike history parameters
V.est_F=0;          % do NOT estimate fluorescence parameters
V.est_c=0;          % do NOT estimate calcium parameters
V.Nspikehist=1;     % assume 1 spike history terms
V.Nparticles=50;    % # of monte carlo samples
V.smc_plot=0;       % do not plot smc
V.smc_iter_max=1;   % number of iterations
V.fast_do=0;        % do not do fast-oopsi
V.smc_do=1;         % do not do smc-oopsi

% set up a few variables
oo = zeros(V.Ncells);
% use default options for glmnet
if ~exist('options','var'), options = glmnetSet; end
% set the number of folds for cross validation
if ~isfield(V,'kFold'), V.kFold = 5; end
kFold = V.kFold;

% estimate spike history terms
hh = zeros(V.Ncells,V.T);
for j=1:V.Ncells
    tau = V.dt / (1-fast{j}.P.gam);
    hh(j,:) = filter(1,[1 -(1-V.dt/tau)],fast{j}.n);
end
% append external input
hh = [V.x; hh];

% estimate connectivity (using netfit)
netinf = cell(V.Ncells,1);
for i=1:V.Ncells
    % get rid of self term
    h = [hh(1:i,:); hh(i+2:end,:)];
    % we're going to cross validate by breaking up the movie into k 
    % sequential chunks of data. thus large k values don't make sense
    folds = 1:length(V.x);
    folds = ceil((folds)/(length(folds)/kFold));
    % argmax w over P(n|h)
    cv  = cvglmnet(h',fast{i}.n,kFold,folds,'response','gaussian',options,0);
    % save best lambda row, determined by k fold cross validation
    cv.lambda_index = find(cv.lambda_min == cv.glmnet_object.lambda);
    % NOTE: BETA contains inferred weights for all lambda
    % however the first index corresponds to the external stimulus term
    % which is all ones by default. the current cell is also missing
    % consequently, to reconstruct a connectivity matrix like we do below
    % you need to add a zero to each row along the diagonal (e.g. if ii is
    % the current row, index ii needs to be added to the middle of the row
    % and set to zero.
    row = cv.glmnet_object.beta(1:end,cv.lambda_index);
    % append the weights to the connectivity matrix
    oo(i,:) = [row(2:i,:); 0; row(i+1:end,:)];
    % save struct
    netinf{i} = cv;
    
end

% save the output
Phat.netinf = netinf;
Phat.omega = oo;
Phat.h = hh;