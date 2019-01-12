dataFileName = fullfile('results', 'doseResponseData.mat');
load(dataFileName);

data = dff;

%% apply PCA, treat each concentration a uniq sample
dataTall = permute(data, [1 3 2]); 
dataTall = reshape(dataTall,[],size(data,2),1);
[coeff, score, ~, ~, explained, mu] = pca (dataTall);

%% plot the data projection onto the first 3 PCs
%define colormap for odors
addpath(genpath(pwd));
[odorColorMap] = FindColorMap();

markerSize = 400;     mksize_temp = [ 0.1 0.2 0.35 0.55 0.8];

figure;
hold on;

for i =1:length(odorList)
	seqTemp = length(odorList) * [0 1 2 3 4] + i;
	scatter3(score(seqTemp, 1), score(seqTemp, 2), ...
        score(seqTemp, 3), markerSize*mksize_temp, odorColorMap(i,:),'fill'); 

	plot3(score(seqTemp, 1), score(seqTemp, 2), ...
        score(seqTemp, 3), 'color', odorColorMap(i,:), ...
        'LineWidth', 1.5);
	pause(0.01);
end

%add odor name
dx = 0.1; dy = 0.1; dz = 0.1; % displacement so the text does not overlay the data points
maxInd = (1:length(odorList)) + 4*length(odorList);
text(score(maxInd, 1) + dx, score(maxInd, 2) +dy, score(maxInd, 3)+dz, odorList);  
xlabel(['PC', num2str(1),' (', num2str(explained(1), 3), '% of variance)']);
ylabel(['PC', num2str(2),' (', num2str(explained(2), 3), '% of variance)']);
zlabel(['PC', num2str(3),' (', num2str(explained(3), 3), '% of variance)']);
axis tight
set(gcf, 'Position', [700 80 1100 850]);
view(159.4, 11.3);
hold off

% save figure
savefig(gcf, fullfile('results', 'figures', 'PCA_Figure2B.fig'));