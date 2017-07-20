load(fullfile('.', 'AnalysisResults', 'FitDoseResponse.mat'));
load(fullfile('.', 'data', 'AveRawDataMatrix.mat'));
data = dataRawAve;
% normalize data of each ORN-odor pair using fitted r_max
[r, c] = find(~isnan(aMatrix));
for i = 1:length(r)
    data(r(i), c(i), :) = data(r(i), c(i), :)/aMatrix(r(i), c(i));
end

data(5, :, :) = [];   odorList(5) = []; %remove odor 5
data(data<0) = 0; %set the negative values to 0

%% define the ONR-to-uPN divisive normalization functions
% Reference: 1.	Olsen, S. R., Bhandawat, V. & Wilson, R. I. 
% Divisive normalization in olfactory population codes. Neuron, 66(2), 287–299 (2010).
rmax = 1;
n = 1.5;
sigma = 0.07;  %12/165
m = 2.5;

normFun = @(x, y) double(rmax .* x.^n./(x.^n + sigma.^n + (m*y).^n));

inputTest = linspace(0, 1, 50);
figure;
co = [ 189,189,189;
    150,150,150;
    115,115,115;
    82,82,82;
    37,37,37];
co = co/255;
LFP = [0.01 0.03 0.06 0.12 0.24];

for i  = 1:length(LFP)
    output= normFun(inputTest, LFP(i));
    plot(inputTest, output, 'Color', co(i,:), 'LineWidth', 1.5 ); hold on;
end
xlabel('ORN'); ylabel('uPN');
title('Transfer Function');
legend(['<ORN>=', num2str(LFP(1))], ['<ORN>=', num2str(LFP(2))], ...
    ['<ORN>=', num2str(LFP(3))], ['<ORN>=', num2str(LFP(4))], ...
    ['<ORN>=', num2str(LFP(5))], 'Location','southeast');

%% normalize and nonlinearize the data
[nOdor, nORN, nConc] = size(data);
PN_Data = zeros(nOdor, nORN, nConc);

for i = 1:nOdor
    for j = 1:nConc
        LFP = mean(data(i, :, j));
        for k = 1:nORN
            PN_Data(i, k, j) = normFun(data(i, k, j) , LFP);
        end
    end
end

%% visualize the ORN data
odorOrder = [17 12 15 2 10 3 4 16 9 1 18 8 6 5 7 13 14 11]; % the order is consistant to figure 2
ORNOrder = [16 17 5 2 14 12 11 1 4 7 10 13 15 9 18 8 6 3];

ORN_DataTemp = data;
[rows, cols, z] = size(data);

ORNnewMStep1 = ORN_DataTemp;
for i = 1:rows
    ORNnewMStep1(i,:,:) = ORN_DataTemp(odorOrder(i),:,:);
end

ORNnewMStep2 = ORNnewMStep1;
for i = 1:cols
    ORNnewMStep2(:,i,:) = ORNnewMStep1(:,ORNOrder(i), :);
end

for i = 1 : z
    data2D = ORNnewMStep2(:,:,i);
    figure;
    a = axes;
    set(a, 'CLim', [0 1]);
    imagesc(data2D); 
    set(gca,'XTick',1:length(infoORNList));
    set(gca,'XTickLabel',infoORNList(ORNOrder));
    set(gca,'xaxisLocation','top');
    set(gca,'YTick',1:length(odorList));
    set(gca,'YTickLabel',odorList(odorOrder));
    ax = gca; ax.XTickLabelRotation = 45;
    title(concList(i));
    colormap(jet);
    caxis([0 1]);
    set(gcf, 'Position', [680   558   510   420]);
end

%% visualize the uPN data
PN_DataTemp = PN_Data;

PNnewMStep1 = PN_DataTemp;
for i = 1:rows
    PNnewMStep1(i,:,:) = PN_DataTemp(odorOrder(i),:,:);
end

PNnewMStep2 = PNnewMStep1;
for i = 1:cols
    PNnewMStep2(:,i,:) = PNnewMStep1(:,ORNOrder(i), :);
end

for i = 1 : z
    data2D = double(PNnewMStep2(:,:,i));
    figure;
    a = axes;
    set(a, 'CLim', [0 1]);
    imagesc(data2D); 
    set(gca,'XTick',1:length(infoORNList));
    set(gca,'XTickLabel',infoORNList(ORNOrder));
    set(gca,'xaxisLocation','top');
    set(gca,'YTick',1:length(odorList));
    set(gca,'YTickLabel',odorList(odorOrder));
    ax = gca; ax.XTickLabelRotation = 45;
    title(concList(i));
    colormap(jet);
    caxis([0 1]);
    set(gcf, 'Position', [680   558   510   420]);
end

%% compare the PCA
% PCA on ORNs
dataTall = permute(data, [1 3 2]); 
dataTall = reshape(dataTall,[],size(data,2),1);
[~,~,~,~,explained_ORN,~] = pca(dataTall );

% plot the variance explained
% figure;
% plot(1:length(explained_ORN), explained_ORN, 'ok');
% xlabel('PC'); ylabel('Variance Explained(%)');
% axis([0 length(explained_ORN) 0 40]);
% title('ORN');

figure;
bar(1:5, explained_ORN(1:5), 0.5, 'k' )
axis([.5 5.5 0 70]);
set(gca,'XTick',1:5);
set(gca,'YTick',0:10:70);
set(gca,'YTickLabel', 0:10:70);
xlabel('Principal Components');
ylabel('% of variance explained');

% PCA on uPNs
dataTallPN = permute(PN_Data, [1 3 2]); 
dataTallPN = reshape(dataTallPN,[],size(PN_Data,2),1);
[~,~,~,~,explainedPN,~] = pca(dataTallPN );

% % plot the variance explained
% figure;
% plot(1:length(explainedPN), explainedPN, 'ok');
% xlabel('PC'); ylabel('Variance Explained(%)');
% axis([0 length(explained_ORN) 0 40]);
% title('uPN');

figure;
bar(1:5, explainedPN(1:5), 0.5, 'k' )
axis([.5 5.5 0 70]);
set(gca,'XTick',1:5);
set(gca,'YTick',0:10:70);
set(gca,'YTickLabel', 0:10:70);
xlabel('Principal Components');
ylabel('% of variance explained');
