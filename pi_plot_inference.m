% plot spike inference(F,I)
% F - cell array of fluorescence vectors
% I - cell array of inference vectors
% s - 1 if neuron is good, 0 otherwise
%
% tamachado 3/10
function pi_plot_inference(F,I,s,t)
figure('Color','w'); hold on;
if ~exist('s','var')
    s = ones(length(F),1);
end
if ~exist('t','var')
    t = sprintf('spike inference (n = %d), excluded in red',length(F));
else
    t = [t sprintf(' (n = %d)',length(F))];
end
rr = 2:2:2*length(F);
for ii = 1:length(F),
    F{ii} = F{ii} - min(F{ii});
    if s(ii), col = 'b'; else col = 'r'; end;
    plot(I{ii} + rr(ii) ,col);
    plot((F{ii}./max(F{ii})) + rr(ii),'k');
end
set(gca,'XTick',[],'XTickLabel',[],...
    'XColor',[1 1 1],'YColor',[1 1 1],...
    'YTick',[],'YTickLabel',[]);
title(t);
orient landscape;
end