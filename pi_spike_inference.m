function [spikes struct] = pi_spike_inference(signals,V,spikes)
% inference = pi_spike_inference(inputs):
%
% inputs (only first argument is required):
%  signal      - fluorescence traces in a matrix (nNeurons x nFrames)
%  V.epath     - path to associated .paq file containing ephys data (path)
%  V.alpha     - threshold for extracting spikes from voltage data (0-1)
%  V.flip      - boolean indicating whether signal should be flipped (boolean)
%  V.channel   - headstage channel containing voltage data (1-4)
%  V.figure    - figure handle to plot stuff in
%  V.indices   - a logical vector showing selected neurons
%  V.fast      - use fast or smc
%  V.inf_plot  - generate spike inference plot
%  spikes      - pass in spike inference here instead of recomputing
%
% output:
%  the choosen spike inference algorithm is run on the input
%  signal. the ephys data is aligned with the fluorescence data.
%  the results are plotted on screen to a new figure.
%
%  inference - vector of length(signal) containing spike inference
%
%  tamachado (5/10)

%% check arguments
% fluorescence data is required
if isempty(signals), error('Fluorescence data must be provided!'); end

% make V if it does not exist
if ~exist('V','var'), V = []; end

% generate plots if not specified
if ~isfield(V,'inf_plot'), V.inf_plot = 1; end

% use fast if not specified
if ~isfield(V,'fast'), V.fast = 1; end

% plotting ephys data is optional
if ~isfield(V,'epath'),   epath = [];      else epath = V.epath; end

% set voltage threshold to default value
if ~isfield(V,'alpha'),   alpha = 0.6;     else alpha = V.alpha; end

% by default, don't flip the traces
if ~isfield(V,'flip'),    flip = false;    else flip = V.flip; end

% add on an arbitrary title
if ~isfield(V,'title'), tt = []; else tt = V.title; end

% did the movie start twice?
if ~isfield(V,'doubleStart'), doubleStart = 0; else doubleStart = V.doubleStart; end

% by default, plot voltage from headstage 1 (with signal from neuron 1)
if ~isfield(V,'channel')
    if exist(epath,'file'), channel = 1; else channel = 0; end
else
    channel = V.channel;
end

% get the number of neurons
Ncells = size(signals,1);
V.Ncells = Ncells;

% update indices structure
if ~isfield(V,'indices')
    indices = ones(Ncells,1);
else
    indices = ones(Ncells,1);
    indices(V.indices) = 0;
end
assignin('base','indices', indices);

% make a new figure if desired
if ~isfield(V,'figure'),  handle = figure; else handle=figure(V.figure); end

%% spike inference parameters

% estimate parameters using smc
V.est_n=1;          % estimate spiking parameters
V.est_h=1;          % estimate spike history parameters
V.est_F=0;          % do NOT estimate fluorescence parameters
V.est_c=0;          % do NOT estimate calcium parameters

% algorithm options
V.fast_iter_max  = 1;                % how many iterations of fast inf
V.fast_plot      = 1;                % whether to generate foopsi plots
V.smc_plot       = 1;                % whether to generate foopsi plots
V.save           = 0;                % whether to save results
V.plot           = 0;                % generate run_oopsi plot
V.Nspikehist     = 1;                % there is 1 history term per neuron
V.StimDim        = 1;                % the stimulus is unidimensional
triggerAmplitude = 3;                % amplitude of camera trigger channel

% needed for generating initial estimates of the spike history term
P.tau_h   = 0.02;                    % decay time constant

if isfield(V,'fast') && V.fast, V.smc_do = 0; V.fast_do = 1; else
                                V.smc_do = 1; V.fast_do = 1; end

%% load up ephys data and imaging fluorescence traces
% find index of camera and voltage signals
cIndex = 0; vIndex = zeros(4,1);
if ~isempty(epath)
    % load up ephys
    info   = paq2lab(epath,'info');
    dataset = paq2lab(epath);
    for ii=1:length(info.ObjInfo.Channel)
        cc = info.ObjInfo.Channel(ii);
        switch cc.ChannelName
            case {'camera','cmera','CameraFrames'}
                cIndex = ii;
            case {'MC1voltage','voltage1','Voltage1'}
                vIndex(1) = ii;
            case {'MC2voltage','voltage2','Voltage2'}
                vIndex(2) = ii;
            case {'MC3voltage','voltage3','Voltage3'}
                vIndex(3) = ii;
            case {'MC4voltage','voltage4','Voltage4'}
                vIndex(4) = ii;
        end
    end
    % give up
    if cIndex == 0, error('Unable to camera index channel!'); end
    if vIndex(1:length(channel) == 0) == 0
        error('Voltage signal on desired channel not found!');
    end
end

% allocate arrays to store spike times
voltage = cell(Ncells,1);
V.ns   = cell(Ncells,1); % spike times from ephys data (in ephys samples)
V.nn   = cell(Ncells,1); % spike times from ephys data (in movie frames)
V.F    = cell(Ncells,1); % fluorescence signals (normalized)
spikes = cell(Ncells,1); % spike inference output vector
struct = cell(Ncells,1); % spike inference output struct
I = cell(Ncells,1);
for cc=1:Ncells,
    I{cc}.M.nbar = zeros(1,size(signals,2));
end

%% process each neuron
for cc = 1:Ncells
    fprintf('\nneuron %d\n',cc);
    % flip the trace if necessary
    if flip == true, fc = -1; else fc =  1; end;
    % extract frame trigger onset times from .paq file
    signal = signals(cc,:);
    % set any zero values to be equal to the mean
    signal(signal == 0) = mean(signal);
    % get the fluorescence trace and flip it if it is fura
    V.F{cc} = fc*signal;
    % preprocess the trace
    V.T    = length(V.F{cc});
    if 1
        V.F{cc} = V.F{cc}-min(V.F{cc}); V.F{cc}=V.F{cc}/std(V.F{cc});
        V.F{cc} = V.F{cc}/max(V.F{cc}); V.F{cc}=V.F{cc}+realmin;
        if ~isfield(V,'flt_opt'), V.flt_opt = 4; end
        switch V.flt_opt
            case 1
                f = V.F{cc};
                % apply filter
                nfft    = 2^nextpow2(V.T);
                y       = fft(f,nfft);
                bw      = 10;
                y(1:bw) = 0; y(end-bw+1:end)=0;
                iy      = ifft(y,nfft);
                V.F{cc}= z1(real(iy(1:V.T)));
                % truncating the end where detrending is bad?
                V.F{cc} = V.F{cc}(1:end-V.flt_rm);
                V.T = length(V.F{cc});
            case 2
                windowSize = 20;
                V.F{cc} =  detrend(medfilt1(V.F{cc},windowSize));
            case 3
                f = detrend(V.F{cc}); ff = f;
                windowSize = 6; delta = .25;
                for ii = 1:length(f)
                    lb = max(ii-windowSize,1); ub = min(ii+windowSize,length(f));
                    mm = median(f(lb:ub));
                    if abs(f(ii)-mm) > delta*mm;
                        ff(ii) = mm;
                    end
                end
                V.F{cc} = ff;
            case 4
                V.F{cc} = detrend(V.F{cc});
        end
    end
    % reload the spike times on each iteration
    times = []; V.nn{cc} = [];
    % what headstage is the corresponding voltage trace on
    [a vC] =  find(channel == cc);
    % load associated ephys
    if ~isempty(vC)
        fprintf('loading associated electrophysiology data...\n');
        % all trigger pulses
        times = find(diff(dataset(:,cIndex)) > triggerAmplitude);
        % display how many there probably are
        fprintf('movie frame times extracted: %d\n',length(times));
        fprintf('actual fluorescence signal length: %d\n',length(signal));
        % electrophysiology sampling rate
        ephysRate = info.ObjInfo.SampleRate;
        % if theres a long off period, extract spikes after it
        if doubleStart
            mTime = find(diff(times) == max(diff(times)));
        end
        % truncate spike time vector to be of fluorescent movie length
        if length(times) > length(signal)
            times = times(1:length(signal));
        end
        % convert from samples to seconds
        times = times * (1/ephysRate);
    end
    % compute frame rate if we have camera traces
    if ~isempty(times)
        V.dt =median(diff(times));
    end
    
    % deal with spike times if desired
    if ~isempty(times)
        % extract spike times from ephys data
        voltage{cc}  = dataset(:,vIndex(vC));
        % if theres a long off period, extract spikes after it
        if doubleStart
            voltage{cc} = voltage{cc}(times(mTime)*ephysRate:end); %#ok<FNDSB>
        end
        V.ns{cc} = get_spike_times(voltage{cc},alpha);
        % convert spike times from ephys samples to movie frames
        nCurr = V.ns{cc};
        for ii = 1:length(nCurr)
            [value nCurr(ii)] = min(abs(times*ephysRate - nCurr(ii)));
        end
        V.nn{cc} = nCurr;
    end
    
    % run the appropriate oopsi
    if isempty(spikes{cc})
        fprintf('running spike inference...\n');
        V.StimDim   = V.Ncells;    % we're fitting coupling terms
        T = V;                     % make a copy of parameters
        T.Ncells = 1;              % there is only 1 neuron in each F
        if V.fast
            struct{cc}    = run_oopsi(V.F{cc},T);
            spikes{cc}    = struct{cc}.n;
            spikes{cc}    = spikes{cc}/max(spikes{cc});
            disp(struct{cc}.P.a/struct{cc}.P.sig);
        else
            % append external stimulus for neuron 'i' with spike histories from other cells
            h = zeros(V.Ncells-1,V.T);          % we append this to x to generate input into neuron from other neurons
            Pre=1:V.Ncells;                     % generate list of presynaptic neurons
            Pre(Pre==cc)=[];                    % remove self
            k=0;                                % counter of dimension
            for dd=Pre                          % loop thru all presynaptic neurons
                k=k+1;                          % generate input to neuron based on posterior mean spike train from neuron j
                h(k,:) = filter(1,[1 -(1-V.dt/P.tau_h)],I{dd}.M.nbar);
            end
            T.x = [V.x; h];
            [out struct{cc}] = run_oopsi(V.F{cc},T);
            % kill useless fields
            struct{cc}.E.w = []; struct{cc}.E.n = []; struct{cc}.E.C = [];
            spikes{cc} = struct{cc}.E.nbar;
        end
    end
end

%% draw the gui
if V.inf_plot
    if ~exist('times','var'), times = 1:length(signal); end;
    figure(handle);
    k=1; kmin=1; kmax=Ncells;
    plot_callback;
    set(gcf,'Color','w','Toolbar','figure');
    guidata(handle,indices);
    if Ncells > 1
    hb = uicontrol(...
        'Style', 'togglebutton',...
        'String', 'Exclude',...
        'Units','normalized',...
        'Position', [0 0 .1 .04],...
        'Callback',@clicked_callback);
    ha = uicontrol(gcf,...
        'Style','slider',...
        'Min' ,kmin,'Max',kmax,...
        'Units','normalized',...
        'Position',[.1 0 .9 .04],...
        'Value', k,...
        'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
        'Callback',@plot_callback);
    % wait until the window is closed before exiting
    uiwait;
    end
end

%% make the interactive plotting window
    function k=plot_callback(handle, eventdata, handles) %#ok<INUSD>
        
        % move the scroll bar
        if exist('ha','var')
            indices = guidata(gcbo);
            k = round(get(ha,'Value'));
        else
            k = 1;
        end
        
        % truncate fluorescence if there are too few camera triggers
        if length(times) < length(V.F{k}) && ~isempty(times)
           V.F{k} = V.F{k}(1:length(times));
           spikes{k} = spikes{k}(1:length(times));
        end
        
        % plot spike times and fluorescence
        ax(1) = subplot(3,1,1);
        cla; plot(V.F{k}./max(V.F{k}),'k');
        mstats = 0;
        if ~isempty(V.nn{k})
            % plot spikes
            hold on; plot(V.nn{k},1,'r.');
            mstats = 1;
            % compute stats
            threshold  = .05;  % significance level
            jitters    = 5;    % area around+- spike that is okay
            espikes = zeros(length(spikes{k}),1);
            espikes(V.nn{k}(V.nn{k} <  length(spikes{k}))) = 1;
            roc = roc3_gamma([spikes{k} espikes],threshold, jitters);
        end
        
        % indicate whether the current neuron is excluded
        if ~exist('indices','var') || length(indices) < k || indices(k) == 1
            tt = sprintf('neuron %d',k);
            eSnr = struct{k}.P.a/struct{k}.P.sig;
            if mstats == 1
                title(sprintf('%s, AUC %.2f, eSNR %.2f',tt,roc.AUC,eSnr),'FontSize',14);
            else
                title(sprintf('%s, eSNR %.2f',tt,eSnr),'FontSize',14);
            end
        else
            title(sprintf('neuron %d: excluded',k),'FontSize',14);
        end
        
        % set x axis appropriately
        nTicks = 40;
        xt = [1:round(length(V.F{k})/nTicks):length(V.F{k})];
        if numel(times) > 0
            xl = round(times(xt) - times(1));
            set(gca,'XTick',xt,'YTick',[],'XTickLabel',xl,'YTickLabel',[]);
            xlabel('time (s)');
        end
        
        % plot inference output
        ax(2) = subplot(3,1,2);
        cla; bar(spikes{k},'k');
        if ~isempty(V.nn{k}), hold on; plot(V.nn{k},1,'r.'); end
        
        
        % link the axes
        linkaxes(ax,'x');
        
        % plot spike times and ephys data if they exist
        subplot(3,1,3);
        if ~isempty(voltage{k}) && ~isempty(V.ns{k})
            cla; plot(voltage{k}./max(voltage{k}),'k');
            hold on; plot(V.ns{k},1,'r.');
            set(gca,'YTick',[],'YTickLabel',[]);
            xlabel('time (s)');
        else
            cla; set(gca,'Visible','off');
        end
    end

    function c=clicked_callback(handle, eventdata, handles) %#ok<INUSD>
        % flip the bit at the appropriate index
        c = get(hb,'Value');
        k = get(ha,'Value');
        % update the data
        indices = guidata(gcbo);
        % save out our modifications
        if indices(k), indices(k) = 0; else indices(k) = 1; end;
        guidata(gcbo,indices);
        plot_callback;
        assignin('base','indices', indices);
    end
end
