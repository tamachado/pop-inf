% run_pi
%
% this script runs population inference on raw fluorescence traces
% and displays the results
%
% in addition to the estimated connectivity matrix, the spike rate
% nonlinearity for each neuron is plotted. this function should be
% increasing (ie the predicted firing rate given the model is correlated
% with the actual firing rate of that neuron). run the code on the 
% simulated data to get an example for what this should look like.
%
% tamachado 5/10

% clear the workspace
clear all;

% set up parameters (or pi_infer_connectivity will prompt for necessary params)
switch 3
    case 1
        % example of how to analyze an arbitrary new dataset
        V.flip = false;
        V.dt = .03333;
        V.path = '../test2.mat';
        V.inf_plot = 1;
    case 2
        % example of how to analyze a dataset that has
        % parameters hardcoded into pi_infer_connectivity
        V.name = 'long-1-of-4';
        V.inf_plot = 1;
        V.all = 1;
    case 3
        % example of how to load in traces from a tif stack compressed
        % using the "compress-stack" function in imagej
        V.flip = true;
        V.dt = .03333;
        V.inf_plot = 1;
        V.traces = pi_load_data([pwd '/data/']);
end

% run connectivity inference
O = pi_infer_connectivity(V);

% plot spike inference summary plot
pi_plot_inference(O.F,O.N,O.indices);

% plot connectivity matrix for lam_max
figure; imagesc(O.Phat.omega); colormap gray;
title('connectivity matrix'); xlabel('postsynaptic'); ylabel('presynaptic');

% plot nonlinearity to show predictive power of model fit
pi_plot_nonlinearity(O);