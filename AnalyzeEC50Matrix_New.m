load('fitResults.mat');

%% load odor solubility
[mw, rho, sol] = loadOdorSolub();
hc = log10(10^-3 * sol./rho);
% cMatrixMol = cMatrix;
% 
% for i = 1:rowNum
%     factor =  10^3 * rho(i) / mw(i);
%     idx = find(cMatrix(i,:) ~= 0);
%     cMatrixMol(i, idx) = cMatrix(i,idx) + log10(factor); 
% end

%% sort the log10(EC50) for each odor
[rowNum, colNum] = size(cMatrix);
cMin = min(cMatrix, [], 2);
[~, rowIdx] = sort(cMin, 'descend');

% show the distribution of the sensitive EC50 for each odor
figure; histogram(cMin);

figure; 
for i = 1:rowNum
    idx = rowIdx(i);
    kdVec = cMatrix(idx, :);
    kdVec = kdVec(find(kdVec));
    [kdVec, ~] = sort(kdVec);   
    
    iShow = i;
    
    yPlot = repmat([iShow-0.3; iShow+0.3], [1, length(kdVec)]);
    xPlot = repmat(kdVec, [2, 1]);
    plot(xPlot, yPlot, 'k'); 
    hold on;
    plot([-10, -2], [iShow, iShow], 'k');
%     plot([hc(idx); hc(idx)], [iShow-0.3; iShow+0.3],'r');
end
xlabel('log_{10}EC_{50}');axis tight; 
% set(gca,'box','off','ycolor','w')
set(gca,'box','off');
set(gcf, 'Position', [200 10 560 988]); title('Each row is an odor');
yticks(1 : rowNum);
yticklabels(odorList(rowIdx));

%% sort the -log10(EC50) for each odor
negC = -cMatrix;
negCMax = max(negC, [], 2);
[~, rowIdx] = sort(negCMax, 'descend');

figure; 
for i = 1:rowNum
    idx = rowIdx(i);
    kdVec = negC(idx, :);
    kdVec = kdVec(find(kdVec));
    [kdVec, ~] = sort(kdVec);   
    
    iShow = i;
    
    yPlot = repmat([iShow-0.3; iShow+0.3], [1, length(kdVec)]);
    xPlot = repmat(kdVec, [2, 1]);
    plot(xPlot, yPlot, 'k'); 
    hold on;
    plot([2, 10], [iShow, iShow], 'k');
%     plot([hc(idx); hc(idx)], [iShow-0.3; iShow+0.3],'r');
end
xlabel('-log_{10}EC_{50}');axis tight; 
% set(gca,'box','off','ycolor','w')
set(gca,'box','off');
set(gcf, 'Position', [200 10 560 988]);
yticks(1 : rowNum);
yticklabels(odorList(rowIdx));


%% sort the EC50 for each ORN
[rowNum, colNum] = size(cMatrix);
cMin = min(cMatrix, [], 1);
[~, colIdx] = sort(cMin);

figure; 
for i = 1:colNum
    idx = colIdx(i);
    kdVec = cMatrix(:, idx);
    kdVec = kdVec(find(kdVec));
    [kdVec, ~] = sort(kdVec');   
    
    iShow = colNum - i + 1;
    
    yPlot = repmat([iShow-0.3; iShow+0.3], [1, length(kdVec)]);
    xPlot = repmat(kdVec, [2, 1]);
    plot(xPlot, yPlot, 'k'); 
    hold on;
    plot([-10, -2], [iShow, iShow], 'k');
end
xlabel('log_{10}EC_{50}'); axis tight; set(gca,'box','off','ycolor','w')
set(gcf, 'Position', [500 10 560 988]); title('Each row is an ORN');

%% plot the distribution of Delta_Kd for each odor
[rowNum, colNum] = size(cMatrix);
cMin = min(cMatrix, [], 2);
[~, rowIdx] = sort(cMin);

figure; 
dKdRowPool = [];
for i = 1:rowNum
    idx = rowIdx(i);
    kdVec = cMatrix(idx, :);    kdVec = kdVec(find(kdVec));
    [kdVec, ~] = sort(kdVec);   
    
    dKd = kdVec(2:end) - kdVec(1);	dKdRowPool = [dKdRowPool dKd];
    iShow = rowNum - i + 1;
    
    yPlot = repmat([iShow-0.3; iShow+0.3], [1, length(dKd)]);
    xPlot = repmat(dKd, [2, 1]);
    plot(xPlot, yPlot, 'k'); 
    hold on;
    plot([0, 7], [iShow, iShow], 'k');
end
xlabel('log_{10}(Kd/Kd_0)');axis tight; set(gca,'box','off','ycolor','w')
set(gcf, 'Position', [200 10 560 988]); title('Each row is an odor');

figure; histogram(dKdRowPool);

