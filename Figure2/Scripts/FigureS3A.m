%%Figure S3A
%PCA analysis and visualization of 32 molecular descriptors for the 690 
%odorants from the DoOr database (gray dots), and the 35 odorants in
%our panel (red dots)

%% clear workspace
clear all;  clc

%% read data
%define and load the .xlsx file
filename = 'eDragon Descriptors for Panel of 35 Odorants';

%Sheet1 of excel file is all 690 odors from DoOr database
allOdorsSheet = 1;
xlRange = 'A1:BLE691'; 
[~,~,raw] = xlsread(filename,allOdorsSheet,xlRange);

%Sheet2 of excel file is the 35 odors from our paper
testOdorsSheet = 2;
x2Range = 'A1:BLE36'; 
[~,~,rawTestOdors] = xlsread(filename,testOdorsSheet,x2Range);

%% Re-format data

%define the information list of the odors
infoOdorAll  = raw(2:end, 1:3);
infoOdorName = infoOdorAll(:, 2);

infoOdorTestOdors  = rawTestOdors(2:end, 1:3);
infoOdorNameTestOdors = infoOdorTestOdors(:, 2);

%define the information list of the descriptors
infoMetric = raw(1, 4:end);

%extract the data matrix
m = cell2mat(raw(2:end, 4:end));

%% PCA on 32 descriptors with axes projections
%list of 32 descriptors as defined in Haddad paper
optList  = [22 48 92 96 359 404	433	513	582	592	678	688	690	698	946 ...
    948	959	963	1012 1069 1110 1191 1286 1321 1331 1348 1373 1430 1528 ...
    1541 1558 1576];
optWeight= [1	1	1	1	1	1	1	1	1	1	1	1	1	1	1 ...
    2	1	3	1	1	1	1	1	1	3	1	2	1	1	2	2	1];

mOpt = m(:, optList);
zOpt = zscore(mOpt); 
[coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(zOpt, 'VariableWeights', optWeight);

%score2 of all odors
x = score2(:, 1);
y = score2(:, 2);
z = score2(:, 3);

%find range of points[min, max, # of points so bin size = 1]
range = [min(x),max(x); min(y),max(y); min(z),max(z)];
range(:,3) = range(:,2) - range(:,1) + 1; %for bin size = 1
%range(:,3) = (range(:,2) - range(:,1) + .5)/.5; %for bin size =0.5

%define the x,y,z bin edges, num is # of points
xi = linspace(range(1,1),range(1,2),range(1,3));
yi = linspace(range(2,1),range(2,2),range(2,3));
zi = linspace(range(3,1),range(3,2),range(3,3));
n = length(xi); m = length(yi); p = length(zi);

%bin points from data set of all odors and normalize
xr = interp1(xi,1:numel(xi),x,'nearest');
yr = interp1(yi,1:numel(yi),y,'nearest');
zr = interp1(zi,1:numel(zi),z,'nearest');
density_odors = accumarray([xr yr zr],1,[n m p]);
odorsAll = density_odors./length(infoOdorName);

%find location of core panel odors (18 odors from Mathews et al, 2013)
panel = zeros(length(score2),1);
loc = cell2mat(infoOdorTestOdors(:,1));
panel(loc,1) = 1;

%score2 for different combinations of panel odors
for i = 1:length(score2)
    if panel(i,1) == 0
        panelNew = panel;
        panelNew(i,1) = 1;
        
        %find score2 of sample odors
        Indx = find(panelNew == 1);
        xSample = score2(Indx, 1);
        ySample = score2(Indx, 2);
        zSample = score2(Indx, 3);
        
        %bin points from new sample odor panel and normalize
        xp = interp1(xi,1:numel(xi),xSample,'nearest');
        yp = interp1(yi,1:numel(yi),ySample,'nearest');
        zp = interp1(zi,1:numel(zi),zSample,'nearest');
        density_panel = accumarray([xp yp zp],1,[n m p]);
        odorsPanel = density_panel./length(xSample);
        
        %find difference matrix between all and sample odors
        odorsDiff = odorsAll - odorsPanel;
        diffMat{i,1} = infoOdorName{i,1};
        diffMat{i,2} = sum(odorsDiff(:));
        
    else
        diffMat{i,1} = infoOdorName{i,1};
        diffMat{i,2} = 0;
    end
    
end

%sort rows according to col 2 = metric value
rankOdors = sortrows(diffMat,2);
    
%% Visualization of 2 PCs with position label for an odor of choice

%pick odor whose position you would like to visualize
odorant = 'trans-3-Hexen-1-ol';
newOdorLoc = find(strcmp(infoOdorName, odorant));

testPanel = panel;
testPanel(newOdorLoc,1) = 1;

%PC1 and PC2
figure; 
h = scatterhist(x,y,'Kernel','on','Group', testPanel,'Color',[80/255,80/255,80/255;255/255,78/255,42/255],...
    'MarkerSize',4.5,'Legend', 'off','LineStyle',{'-','-'});
h(1).Children(1).MarkerEdgeColor = 'black';
h(1).Children(1).MarkerFaceColor = [255/255,78/255,42/255];
h(1).Children(2).MarkerEdgeColor = [99/255,99/255,99/255];
h(1).Children(2).MarkerFaceColor = [220/255,220/255,220/255];

title('PC1 and PC2');
xlabel(['First principal component (', num2str(explained2(1), 3), '% of variance)']);
ylabel(['Second principal component (', num2str(explained2(2), 3), '% of variance)']);
%text for odor
 text(score2(newOdorLoc,1) +.1, score2(newOdorLoc,2) +.3, infoOdorName(newOdorLoc,1));
% % text for all DoOR odors
%   text(x+.1, y+.2, infoOdorName);

%PC1 and PC3
figure; 
h = scatterhist(x,z,'Kernel','on','Group', testPanel,'Color',[80/255,80/255,80/255;255/255,78/255,42/255],...
    'MarkerSize',4.5,'Legend', 'off','LineStyle',{'-','-'});
h(1).Children(1).MarkerEdgeColor = 'black';
h(1).Children(1).MarkerFaceColor = [255/255,78/255,42/255];
h(1).Children(2).MarkerEdgeColor = [99/255,99/255,99/255];
h(1).Children(2).MarkerFaceColor = [220/255,220/255,220/255];

title('PC1 and PC3');
xlabel(['First principal component (', num2str(explained2(1), 3), '% of variance)']);
ylabel(['Third principal component (', num2str(explained2(3), 3), '% of variance)']);
%text for new odor
 text(score2(newOdorLoc,1) +.1, score2(newOdorLoc,3) +.3, infoOdorName(newOdorLoc,1));
% % text for all odors
% text(x+.1, z+.1, infoOdorName);

%PC2 and PC3
figure; 
h = scatterhist(y,z,'Kernel','on','Group', testPanel,'Color',[80/255,80/255,80/255;255/255,78/255,42/255],...
    'MarkerSize',4.5,'Legend', 'off','LineStyle',{'-','-'});
h(1).Children(1).MarkerEdgeColor = 'black';
h(1).Children(1).MarkerFaceColor = [255/255,78/255,42/255];
h(1).Children(2).MarkerEdgeColor = [99/255,99/255,99/255];
h(1).Children(2).MarkerFaceColor = [220/255,220/255,220/255];

title('PC2 and PC3');
xlabel(['Second principal component (', num2str(explained2(2), 3), '% of variance)']);
ylabel(['Third principal component (', num2str(explained2(3), 3), '% of variance)']);
%text for new odor
text(score2(newOdorLoc,2) +.1, score2(newOdorLoc,3) +.3, infoOdorName(newOdorLoc,1));
% % text for all odors
% text(y+.1, z+.3, infoOdorName);

%% Visualization of 3 PCs (Fig. S3A)
figure; 
S = testPanel*8 + 8;
C = -1*(testPanel-1).* ones(690, 3) * 99/255 + ...
    testPanel .* [testPanel, repmat(78/255, 690, 1), repmat(42/255, 690, 1)] ;
scatter3(x,y,z,S,C, 'filled');
xlabel('PC1');  ylabel('PC2');  zlabel('PC3');
