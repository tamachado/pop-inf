function pi_plot_nonlinearity(O)
%% set stuff up
Ncells = size(O.Phat.omega,1);
figure;
k=1; kmin=1; kmax=Ncells;
plot_callback;
set(gcf,'Color','w','Toolbar','figure');
if Ncells > 1
    ha = uicontrol(gcf,...
        'Style','slider',...
        'Min' ,kmin,'Max',kmax,...
        'Units','normalized',...
        'Position',[.1 0 .9 .04],...
        'Value', k,...
        'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
        'Callback',@plot_callback);
end

%% make the interactive plotting window
    function k=plot_callback(handle, eventdata, handles)
        % move the scroll bar
        if exist('ha','var')
            ii = round(get(ha,'Value'));
        else
            ii = 1;
        end
        
        % get generator function (predicted spike inference)
        object = O.Phat.netinf{ii}.glmnet_object;
        nbeta=[object.a0'; object.beta];
        nbeta=nbeta(:,O.Phat.netinf{ii}.lambda_index);
        % since self terms were removed from h during fitting, we must add
        % a zero in at the appropriate location to make the indexing work
        % we add the zero at ii+1 because we added a0 to nbeta already,
        % where a0 is the intercept
        nbeta =  [nbeta(1:ii+1,:); 0; nbeta(ii+2:end,:)];
        switch 2
            case 1
                % predicted spike inference
                k = [ones(size(O.Phat.h',1),1), O.Phat.h'] * nbeta;
            case 2
                % predicted spike inference (exclude current cell)
                k = [ones(size(O.Phat.h',1),1), O.Phat.h(1:ii,:)',O.Phat.h(ii+2:end,:)']...
                    * [nbeta(1:ii); nbeta(ii+2:end)];
        end
        
        % get spike inference (convolved with exponential)
        hc = O.Phat.h(ii+1,:);
        
        % how much variance is explained by the model
        %hcn = hc./max(hc);
        res = k - hc';
        vv = (res'*res)/(hc*hc');

        
        % plot nonlinearity
        subplot(1,2,1); cla; title(sprintf('Neuron %d, residual variance = %.2f', ii, vv));
        plot(k,hc,'.','Color',[0.5725 0.5725 0.5725])
        nBins = 100;
        x = zeros(nBins-1,1);
        y = x; e = x;
        percent = linspace(0,100,nBins);
        for jj = 1:length(percent)-1;
            ind = find(k >= prctile(k,percent(jj)) & ...
                k < prctile(k,percent(jj+1)));
            y(jj) = prctile(hc,mean([percent(jj) percent(jj+1)]));
            x(jj) = mean(k(ind));
            e(jj) = std(k(ind));
        end
        hold on; errorbar(x,y,e,'k','LineWidth',3);
        xlabel('generator signal');
        ylabel('spike inference');
        
        % plot spike trains
        subplot(1,2,2); cla;
        plot(k./max(k) - 1,'b'); 
        hold on;
        plot(hc ./ max(hc) + 1,'r');
        plot((k./max(k)) - (hc ./ max(hc))','k');
        legend('generator','data','residual');
        set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);
    end
end
