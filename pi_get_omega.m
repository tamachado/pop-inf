%
% pass in output from pr_infer_connectvity_condor
% get omega at the specified value for lambda (min, max, best)
% this function assumes there is one external stimulus variable
%
function omega = pi_get_omega(out,lam,off)
% change this if there is more than 1 external stimulus variable
nStimExt = 1;
% allocate space for the connectivity matrix
omega = zeros(out.V.Ncells);
% offset from min/max lambda
if ~exist('off','var'), off = 1; end
for i=1:out.V.Ncells
    beta = out.Phat.netinf{i}.beta;
    % get specified row
    switch lam
        case 'min'
            row = -beta(1+nStimExt:end,1+off);
        case 'max'
            row = -beta(1+nStimExt:end,end-off);
        case 'mean'
            row = -mean(beta(1+nStimExt:end,:),2);
        case 'median'
            row = -mean(beta(1+nStimExt:end,:),2);
        case 'best'
            error('no');
    end
    % append weights to growing matrix
    omega(i,:) = row;
end