% addpath(genpath(pwd));
% [dataMean, dataSEM, odorList, concList, infoORNList ] = GetAveNormData(); %get the averaged and normalized data
% 
% data = dataMean; 
% dataSEM(5,:,:)=[]; 

% load(fullfile('.', 'AnalysisResults', 'FitDoseResponse.mat'));

% load(fullfile('.', 'data', 'AveRawDataMatrix.mat'));
dataFileName = fullfile('.', 'data', 'AveRawDataMatrix2ndRound.mat');
load(dataFileName);

data = dataMean;

% % normalize data of each ORN-odor pair using fitted r_max
% [r, c] = find(~isnan(aMatrix));
% for i = 1:length(r)
%     data(r(i), c(i), :) = data(r(i), c(i), :)/aMatrix(r(i), c(i));
% end

% data(5, :, :) = [];   odorList(5) = []; %remove odor 5

%% apply PCA, treat each concentration a uniq sample
dataTall = permute(data, [1 3 2]); 
dataTall = reshape(dataTall,[],size(data,2),1);
[coeff, score, ~, ~, explained, mu] = pca (dataTall);

center = mean(dataTall, 1);

%% plot the variance explained
figure; 
plot(1:length(explained), explained, 'ok');
axis([0 length(explained) 0 40]);
xlabel('PC'); ylabel('Variance Explained(%)');

% axis([0 length(explained) 0 30]);

%% plot the data projection onto the first 3 PCs
%define colormap for odors

[odorColorMap] = FindColorMap();

% odorColorMap = jet(34);

% odorColorMapRaw = [...
%     138 198 64;...	%1-pentanol, 35a
%     175 52 147 ;...	%3-pentanol
%     0 145 87 ;...	%6-methyl-5-hepten-2-ol, 13a
%     78 71 157;...	%3-octanol, 85c
%     251 106 74;...	%methyl phenyl sulfide, 24a
%     239 59 44; ...	%anisole, 30a
%     203 24 29 ; ...	%2-acetylpyridine, 22c
%     217,95,14; ...  %2,5-dimethylpyrazine, 33b
%     67 182 73; ...  %pentyl acetate, 47a
%     118,42,131  ;...%geranyl acetate, 82a
%     140,81,10 ; ... %2-methoxyphenyl acetate, 94a
%     45 85 166; ...  %trans,trans-2,4-nonadienal, 
%     103 0 13;...    %4-methyl-5-vinylthiazole, 45b
%     165 15 21; ...  %4,5-dimethylthiazole, 59a
%     64,0,75; ...    %4-hexen-3-one, 42a
%     30 142 205; ...	%2-nonanone, 45a
%     153,112,171; ...%acetal, 42b
%     254 192 15; ... %2-phenyl ethanol
%     ];
% odorColorMap = odorColorMapRaw./255.0;
% odorColorMap = zeros(length(odorList), 3);
% odorColorMap(1:19 , :) = repmat([0 0 1], [19, 1]);
% odorColorMap(20:length(odorList) , :) = repmat([1 0 0], [length(odorList) - 19, 1]);
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
view([1,2,1]);
hold off

% %% calculate the distance and angle in the 3D space
% concenN = 5; %number of concentrations
% 
% % formate the data to 3D matrix
% dataPCA = score(:, 1:3);
% [rows, cols] = size(dataPCA);
% newRows = rows/concenN;
% 
% dataPCA3D = zeros(newRows, cols, concenN);
% for i = 1 : newRows
%     for j = 1:concenN
%         index = (j-1)*newRows + i;
%         dataPCA3D(i, :, j) = dataPCA(index, :);
%     end
% end
% 
% % define the original point as the center of the lowest concentration data
% origin = -center;
% originPCA = origin*coeff; 
% originPCAM = repmat(originPCA(1:3), newRows, 1, concenN);
% 
% %shift the center to the origin
% dataRelative = dataPCA3D - originPCAM;
% 
% % calculation to 
% startConc = 1;
% r_theta_phi = zeros(newRows, cols, concenN);
% for i = 1:newRows
%     for j = startConc:concenN
%         x = dataRelative(i,1,j);
%         y = dataRelative(i,2,j);
%         z = dataRelative(i,3,j);
%         r_theta_phi(i,1,j) = sqrt(x^2 + y^2 + z^2); %radius, r
%         r_theta_phi(i,2,j) = acos(z/sqrt(x^2 + y^2 + z^2)); %inclination, theta
%         r_theta_phi(i,3,j) = atan(y/x); %azimuth, phi
%     end
% end
% 
% %% plot
% % plot phi vs theta
% figure;
% markerSize2= 5;
% mksizeFactor2 = [ 0.2 0.4 0.8 1.1 1.4];
% for i = 1:newRows
%     for j = 3:concenN
%         plot(r_theta_phi(i, 2, j)/pi,r_theta_phi(i, 3, j)/pi, ...
%             'Color', odorColorMap(i,:), ...
%             'Marker', 'o',...
%             'Color',  odorColorMap(i,:), ...
%             'MarkerSize', markerSize2*mksizeFactor2(j), ...
%             'MarkerEdgeColor', odorColorMap(i,:), ...
%             'MarkerFaceColor', odorColorMap(i,:) );
%         hold on;
%     end
% end
% hold off;
% axis([0 1 -0.5 0.5]);
% ax = gca; 
% ax.XTick = [0  0.25  0.5  0.75  1]; 
% ax.YTick = [-0.5 -0.25 0 0.25 0.5];
% ax.XTickLabel = {'0', '1/4', '1/2', '3/4', '1'}; 
% ax.YTickLabel = {'-1/2', '-1/4', '0', '1/4', '1/2'};
% xlabel('\Theta');
% ylabel('\Phi');
% 
% %plot the r vs concentration
% figure;  
% inten = zeros(1, 5);
% for i =1:length(odorList)
%     for j = 1:length(concList)
%         inten(j) =  r_theta_phi(i, 1, j);
%     end
%     plot([1:5] + 0.02*i, inten, 'o', 'MarkerEdgeColor', odorColorMap(i, :), 'MarkerFaceColor', ...
%         odorColorMap(i, :), 'MarkerSize',6); 
%     hold on;
% end
% hold off; axis([1  5.5  -0.1  2]); 
% ax = gca; 
% ax.XTick = [1.25  2.25  3.25  4.25  5.25]; 
% ax.YTick = [0 0.5 1 1.5 2];
% ax.XTickLabel = {'1E-8', '1E-7', '1E-6', '1E-5', '1E-4'}; 
% ax.YTickLabel = {'0', '0.5', '1', '1.5', '2'};
% xlabel('Concentrations'); ylabel('ORNs total activities (a.u.)');