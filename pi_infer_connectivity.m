% population_inference = pi_infer_connectivity(S)
%
% infer spikes using fast-oopsi and then
% infer connectivity on these fluorescence traces
%
% arguments:
% S is a struct with the following fields
% either S.name or S.path must be provided
% S -- name: name of dataset to analyze
%   -- path: path to dataset to use
%
% the dataset file must contain a matrix in the form nCells x nFrames
% this matrix must be named "traces"
%
% look at the default parameters set inside the function
% to see what other parameters can be specified
function O = pi_infer_connectivity(V)
try
    tic
    %% choose a dataset
    
    %%%%%%%%%%%%%%%%%%%%
    % error checking
    %%%%%%%%%%%%%%%%%%%%
    if ~exist('V','var'), error('A struct must be passed in'); end
    if ~isfield(V,'name'), V.name = 'null'; end
    % always use fast spike inference
    V.fast = true;
    pf = 'fast';
    % don't draw spike inference gui by default
    if ~isfield(V,'inf_plot'), V.inf_plot = 0; end
    traces = [];
    
    %%%%%%%%%%%%%%%%%%%%
    % default parameters
    %%%%%%%%%%%%%%%%%%%%
    if ~isfield(V,'save'),    V.save   = 0;  end     % save output mat file
    if ~isfield(V,'kFold'),   V.kFold = 5;   end     % sequential chunks to make for cross validation
    if ~isfield(V,'StimDim'), V.StimDim = 1; end     % dimensionality of stimulus
    if ~isfield(V,'indices'), V.indices = [];end     % excluded neurons
    if ~isfield(V,'lam'),     V.lam     = 0; end     % expected firing rate (see fast-oopsi)
    if ~isfield(V,'est_sig'), V.est_sig = 0; end     % sd of noise (see fast-oopsi)
    if ~isfield(V,'est_gam'), V.est_gam = 0; end     % calcium decay tau (see fast-oopsi)
    if ~isfield(V,'est_lam'), V.est_lam = 0; end     % expected firing rate (see fast-oopsi)
    if ~isfield(V,'est_a'),   V.est_a   = 1; end     % spatial filter (see fast-oopsi)
    if ~isfield(V,'est_b'),   V.est_b   = 1; end     % background (see fast-oopsi)          
    options   = glmnetSet;                           % use default glmnet options
    if ~isfield(V,'flt_opt'), V.flt_opt = 4; end     % which filter to use on traces:
    % look at pi_spike_inference for more details
                    % 1 = bandpass using inverse fft
                    % 2 = detrend + matlab median filter
                    % 3 = detrend + rolling median filter
                    % 4 = detrend only
    if ~isfield(V,'flt_rm'), V.flt_rm = 500; end
    % if we filter using the bandpass filter, we remove this many samples
    % at the end because of weird aliasing effects.
    
    %%%%%%%%%%%%%%%%%%%%
    % choose the dataset
    %%%%%%%%%%%%%%%%%%%%
    switch V.name
        % which quarter of the data should be used
        case 'long-1-of-4'
            p = 1;
        case 'long-2-of-4'
            p = 2;
        case 'long-3-of-4'
            p = 3;
        case 'long-4-of-4'
            p = 4;
        case 'long-all'
            p = 0;
        case 'test'
            load('../test2.mat');
            %traces = traces(1:10,:);
            V.dt      = 0.0333;       % 1/frame_rate
            V.flip    = false;        % flip fluorescence traces
        case 'slm'
            load('../slm.mat');
            V.dt      = 0.01567;      % 1/frame_rate
            V.flip    = true;         % flip fluorescence traces
        case 'slm2'
            load('../slm2.mat');
            V.dt      = 0.01567;      % 1/frame_rate
            V.flip    = true;         % flip fluorescence traces
            V.flt_opt = 1;            % bandpass filter
        otherwise
            if isfield(V,'traces')
                traces = V.traces;
            else
                if ~isfield(V,'path')
                    V.path = input('Path to data: ', 's');
                end
                load(V.path);
            end
            if ~isfield(V,'dt')
                V.dt = input('Input frame length in seconds: ');
            end
            if ~isfield(V,'flip')
                V.flip = input('Flip fluorescence traces (1 or 0): ');
            end     
    end
    % create the external stimulus if it hasn't been loaded elsewhere in
    % order to normalize for the effect of differental mean firing rates
    if ~isfield(V,'x'), V.x = ones(V.StimDim,length(traces)); end
    % if we're using the long in vivo (vb-100117) dataset, load it up
    if findstr('long',V.name)
        load('../timecourses2.mat');
        traces = timeCourses';
        clear timeCourses;
        % use only a quarter of the data
        if p > 0
            q = length(traces)/4;
            traces = traces(:,(p-1)*q+1:p*q);
        end
        % use all neurons or just good neurons
        if ~isfield(V,'all') || ~V.all
            V.indices = ...
                [1 2 4 7 8 9 13 16 23 35 36 38 39 47 51 55 64 67 70 78 79];
        else
            V.indices = [];
        end
        % set glm options appropriately
        options.standardize = 1;
        options             = glmnetSet(options);
        % set other options
        V.dt        = 0.0326;
        V.flip      = false;
        V.StimDim   = 1;
        V.x         = ones(V.StimDim,length(traces));
        V.lam       = 10;
    end

    %% run spike inference
    [spikes S] = pi_spike_inference(traces,V);
    indices = ones(size(traces,1),1);
    indices(V.indices) = 0;

    %% update variables appropriately
    V.name = [V.name pf];
    validIndices = find(indices == 1);
    F = cell(length(validIndices),1);
    SS = cell(length(validIndices),1);
    V.Ncells = length(F); V.T=size(traces,2); V.StimDim = V.Ncells;
    if V.flt_opt == 1
        V.T = V.T - V.flt_rm;
        V.x = V.x(1:V.T);
    end
    for ii = 1:V.Ncells
        F{ii}  = traces(validIndices(ii),:);
        SS{ii}.P = S{validIndices(ii)}.P;
        SS{ii}.n = S{validIndices(ii)}.n;
        S{validIndices(ii)} = [];
    end

    % normalize fast spike inference s.t. the largest value is 1
    for ii = 1:length(SS)
        SS{ii}.n = SS{ii}.n ./ max(SS{ii}.n);
    end

    %% use inferred spike trains to infer connectivity
    % get the connectivity matrix
    Phat = pi_run_inference_fast(SS,F,V,options);
    toc
    whos
    
    %% save out stuff
    for ii = 1:length(F), O.N{ii} = SS{ii}.n; end; O.Phat = Phat;
    O.F = F; O.V = V; O.indices = indices;
    if V.save, save([V.name '-' pf],'O'); end
    
catch db
    disp('connectivity inference failed! dumping debug struct (O.db)...')
    whos
    toc
    O.db = db;
end

