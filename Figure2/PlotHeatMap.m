load(fullfile( 'results', 'doseResponseData.mat'));

reproduceFig2 = 1;
% Simulated annealing algorithm is used to search a combination of odor and
% ORN order to best visaulize the heatmap.

if reproduceFig2
    % Because SA is a stochastic algorithm, the result may be different for each 
    % run. To reproduce the figure in the paper, use the following orders 
    % (which gives a loss of 9.54):
    ORNOrder = [19,21,3,6,8,14,10,16,9,7,11,4,12,13,1,2,20,5,17,15,18];
    odorOrder = [19,33,12,32,27,15,14,7,8,9,6,22,24,31,30,29,1,5,10,3,23,4,20,28,26,17,13,25,16,2,21,11,34,18];
else
    % Use simulated annealing (SA) algorithm to search a order of row and col
    [odorOrder, ORNOrder] = GetBestOrder( dffHm );
end

% organize the matrix using the odor and ORN order
mNewOrderTemp = dffHm;
for i = 1:length(odorOrder)
    mNewOrderTemp(i,:,:) = dffHm(odorOrder(i),:,:);
end
mNewOrder = mNewOrderTemp;
for i = 1:length(ORNOrder)
    mNewOrder(:,i,:) = mNewOrderTemp(:,ORNOrder(i), :);
end

%% draw the heatmap
figure; set(gcf, 'Position', [100 300 2870 700] );
for i = 1 : length(concHm)
    data2D = mNewOrder(:,:,i);
    
    a = subplot(1, length(concHm), i);
    
    set(a, 'FontName', 'Arial');
    set(a, 'CLim', [min(min(min(mNewOrder))) 1]);
    imagesc(data2D); 
    set(gca,'XTick',1:length(ORNOrder));
    set(gca,'XTickLabel',ORNList(ORNOrder));
    set(gca,'xaxisLocation','top');
    
    ax = gca; ax.XTickLabelRotation = 90;
    
    if i == 1
        set(gca,'YTick',1:length(odorOrder));
        set(gca,'YTickLabel',odorList(odorOrder));
    else
        set(gca,'ytick',[]);
    end
    
    title(num2str(concHm(i)));
    colormap(jet);
    caxis([min(min(min(mNewOrder))) max(max(max(mNewOrder)))]);

end

% save figure
savefig(gcf, fullfile('results', 'figures', 'heatmap_Figure2A.fig'));

%% plot alcohol group heatmap, Figure S3C
alcoholORN = [4, 13, 11, 12];
alcoholodor= [1,  3,  5,  4];

concNew = concHm(2:end);
alcoholData = dffHm(alcoholodor, alcoholORN, 2:end);

figure;
for i = 1:4
    data = squeeze(alcoholData(i, :, :));
    
    a = subplot(1, 4, i);
    
    set(a, 'FontName', 'Arial');
    set(a, 'CLim', [min(alcoholData(:)) max(alcoholData(:))]);
    imagesc(data); 
    
    set(gca,'XTick',1:4);
    set(gca,'XTickLabel', num2str(concNew));
    
    set(gca,'YTick',1:4);
    set(gca,'YTickLabel', ORNList(alcoholORN));

    
    title(odorList(alcoholodor(i)));
    colormap(jet);
    caxis([min(min(min(alcoholData))) max(max(max(alcoholData)))]);
end
set(gcf, 'Position', [50 500 1650 300]);

% save the figure
savefig(gcf, fullfile('results', 'figures', 'heatmap_FigureS3C.fig'));