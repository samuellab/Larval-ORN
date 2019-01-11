% this code aims to understand what is the main factor causing ORN
% population activity in PCA space roughly a vector

%% setup 
% fix EC50 matrix, tune the slope of the activation function
load('\\LABNAS100\Guangwei\code\fitting\FinalFittingResult\MLEFit_Final_20180724.mat', 'cMatrixMLE');
load('cleanedData.mat', 'odorList');

% replace NaN with 0
cM = cMatrixMLE; cM(isnan(cM)) = 0;

% define x and function form
x = [-8 -7 -6 -5 -4]; 
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

%% case 1, ignore the variations of the 
% uniform amplitude and fixed slop as measured
% show different 
AM = ones(size(cMatrixMLE));
b = 1.56;
y = zeros([size(cM), length(x)]);

for i = 1:size(cM, 1)
    for j = 1:size(cM, 2)
        y(i, j, :) = hillEq(AM(i, j), b, cM(i, j), x); 
    end
end

% PCA
fPlotPCA(y, odorList);

%% case 2, change the slop to a fix but different value
bVec = [0.7 3];
y = zeros([size(cM), length(x)]);

for k = 1 : length(bVec)
    b = bVec(k);
    for i = 1:size(cM, 1)
        for j = 1:size(cM, 2)
            y(i, j, :) = hillEq(AM(i, j), b, cM(i, j), x); 
        end
    end
    
    % PCA and plot
    fPlotPCA(y, odorList);
end

%% case 3, varying slope to be different values among ORNs
bM = rand(size(AM)) * 3.5 + 0.5; %generate uniform distribution between [0.5 4].

y = zeros([size(cM), length(x)]);

for i = 1:size(cM, 1)
    for j = 1:size(cM, 2)
        b = bM(i, j);
        y(i, j, :) = hillEq(AM(i, j), b, cM(i, j), x); 
    end
end

% PCA and plot
fPlotPCA(y, odorList);

%% case 4, shuffle the c-Matrix, remove the correlations
b = 1.56;
cM_Shuffle = zeros(size(cM));

for i = 1:size(cM,1)
    %random shuffle the elements in each row
    temp = cM(i, :);
    cM_Shuffle(i, :) = temp(randperm(size(cM, 2)));
end

y = zeros([size(cM), length(x)]);

for i = 1:size(cM, 1)
    for j = 1:size(cM, 2)
        y(i, j, :) = hillEq(AM(i, j), b, cM_Shuffle(i, j), x); 
    end
end

% PCA and plot
fPlotPCA(y, odorList);

%% case 5_1, change the distribution of the values and check if shape change
% option 1, keep the sparseness and rank, change the distribution to
% uniform distribution 

cM_exp = zeros(size(cM));

cMax = max(cMatrixMLE(:));  cMin = -8;
% cMin = min(cMatrixMLE(:));

for i = 1:size(cM,1)
    %random shuffle the elements in each row
    temp = cMatrixMLE(i, :);
    
    %find out the non-NaN elements
    index = find(~isnan(temp));
    [B, I] = sort(temp(index));
    
    % generate unifrom-distributed data with the boundary 
    newData = rand(size(B)) * (cMax - cMin) + cMin;    
    temp(index(I)) = sort(newData);
    
    cM_exp(i, :) = temp;
end

cM_exp(isnan(cM_exp)) = 0;

y = zeros([size(cM), length(x)]);

for i = 1:size(cM, 1)
    for j = 1:size(cM, 2)
        y(i, j, :) = hillEq(AM(i, j), b, cM_exp(i, j), x); 
    end
end

% PCA and plot
fPlotPCA(y, odorList);

%% case 5_2, change the distribution of the values and check if the shape change
% option 2, keep the known rank, remove sparseness constrain, change the 
% distribution to uniform distribution 

cM_exp = zeros(size(cM));

cMax =0; cMin = min(cMatrixMLE(:));

for i = 1:size(cM,1)
    %random shuffle the elements in each row
    temp = cMatrixMLE(i, :);
    
    %find out the non-NaN elements
    index = find(~isnan(temp));
    [B, I] = sort(temp(index));
    
    % generate unifrom-distributed data with the boundary 
    newData = rand(21,1) * (cMax - cMin) + cMin;    
    newDataTop = mink(newData, length(I));
    temp(index(I)) = sort(newDataTop);
    
    dataCommon = intersect(newData, newDataTop);
    dataLeft = setxor(newData, dataCommon);
    index2 = find(isnan(temp));
    temp(index2) = dataLeft(randperm(length(dataLeft)));
    
    cM_exp(i, :) = temp;
end

cM_exp(isnan(cM_exp)) = 0;

y = zeros([size(cM), length(x)]);

for i = 1:size(cM, 1)
    for j = 1:size(cM, 2)
        y(i, j, :) = hillEq(AM(i, j), b, cM_exp(i, j), x); 
    end
end

% PCA and plot
fPlotPCA(y, odorList);


