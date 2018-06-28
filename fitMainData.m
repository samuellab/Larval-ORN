fileName = fullfile('data', 'Supplementary Table 3.csv');

% load Excel file
dataT = readtable(fileName);

[input.Odor, ~, input.indOdor] = unique(dataT.Odor, 'stable');	%ORN name
input.ORN = dataT.Properties.VariableNames(4:end);
input.concList = dataT.Concentration;
input.expID = dataT.Exp_ID;
input.dff = table2array(dataT(:, 4:end));

%% plot the response curve of every odor 
plotFlag = 0;
if plotFlag == 1
for i = 1 : length(input.Odor)
    figure('Name', input.Odor{i}, 'NumberTitle','off', ...
        'units','normalized','outerposition',[0 0 1 1]);
    hold on;
    
    odI = find(input.indOdor == i); % row index of each odor
    
    concList = input.concList(odI); % concentration list
    dffPool = input.dff(odI, :);    % response data block
    
    [expList, ia, ic] = unique(input.expID(odI), 'stable');   % find out how many trials each odor
    
    for j = 1: length(expList)
        expInd = find(ic == j);
        c = concList(expInd);
        r = dffPool(expInd, :);

        for k = 1:length(input.ORN)
            subplot(3, 7, k); 
            plot(log10(c), r(:, k), 'o-'); title(input.ORN{k});
            xlabel('log10(c)'); ylabel('\DeltaF/F'); hold on;
            pause(0.0001); 
        end
    end
    hold off;
    FigName   = [num2str(get(gcf, 'name')), '.fig'];
    savefig(gcf, fullfile('AnalysisResults', FigName)); close(gcf);
end
end

%% manually set mark good data
maskPlot = zeros(length(input.Odor), length(input.ORN));

maskPlot(1, 4) = 1;     %1-pentanol 35a; S
maskPlot(2, 5) = 1;     %3-pentanol 42a; ?
maskPlot(3, 12) = 1;    %6-methyl-5-hepten-2-ol 85c; S
maskPlot(3, 13) = 1;    %6-methyl-5-hepten-2-ol 13a; S
maskPlot(4, 1)  = 2;    %3-octanol 33b; 
maskPlot(4, 12) = 1;    %3-octanol 85c; S
maskPlot(4, 13) = 1;    %3-octanol 13a; S
maskPlot(5, 4) = 1;     %trans-3-hexen-1-ol 35a; S
maskPlot(6, 8) = 1;     %methyl phenyl sulfide 45b; S
maskPlot(6, 10) = 1;    %methyl phenyl sulfide 24a; S
maskPlot(6, 14) = 2;    %methyl phenyl sulfide 30a;
maskPlot(6, 16) = 1;    %methyl phenyl sulfide 22c; S
maskPlot(8, 10) = 1;    %2-acetylpyridine 24a; S
maskPlot(10, 1) = 1;    %pentyl acetate 33b; S
maskPlot(10, 4) = 1;    %pentyl acetate 35a; S
maskPlot(10, 12) = 1;   %pentyl acetate 85c; S
maskPlot(10, 13) = 1;   %pentyl acetate 13a;
maskPlot(11, 2) = 1;    %geranyl acetate 45a; S
maskPlot(11, 15) = 1;   %geranyl acetate 82a; S
maskPlot(13, 2) = 2;    %trans,trans-2,4-nonadienal 45a; 
maskPlot(13, 20) = 1;   %trans,trans-2,4-nonadienal 74a; S
maskPlot(14, 6) = 1;    %4m5v 59a; S
maskPlot(14, 8) = 1;    %4m5v 45b; S
maskPlot(15, 6) = 2;    %4,5-dimethylthiazole 59a; not flat
maskPlot(16, 5) = 1;    %4-hexen-3-one 42a; S
maskPlot(20, 12) = 1;	%butyl acetate 85c; S
maskPlot(21, 17) = 1;	%ethyl acetate 42b; S
maskPlot(22, 8) = 1;	%benzaldehyde 45b; S
maskPlot(22, 10) = 1;	%benzaldehyde 24a; S
maskPlot(23, 1) = 1;	%2-heptanone 33b/47a; S
maskPlot(23, 4) = 1;	%2-heptanone 35a;  S
maskPlot(23, 5) = 1;	%2-heptanone 42a;  S
maskPlot(23, 12) = 1;	%2-heptanone 85c; 
maskPlot(23, 13) = 2;	%2-heptanone 13a; 
maskPlot(24, 9) = 2;	%methyl salicylate 63a;
maskPlot(24, 10) = 1;	%methyl salicylate 24a; S
maskPlot(24, 12) = 2;	%methyl salicylate 85c; % have strange low point
maskPlot(24, 16) = 1;	%methyl salicylate 22c; S
maskPlot(26, 1) = 1;	%isoamyl acetate 33b/47a; S
maskPlot(26, 12) = 1;	%isoamyl acetate 85c; S
maskPlot(26, 16) = 2;	%isoamyl acetate 22c; % have strange low point
maskPlot(28, 1) = 1;	%hexyl acetate 33b/47a;  S
maskPlot(28, 2) = 2;	%hexyl acetate 45a;  S
maskPlot(28, 4) = 1;	%hexyl acetate 35a;  S
maskPlot(34, 4) = 1;	%nonane 35a; S


%% setup fitting function and method
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));
ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [-1 -1 -11]; % setup the range and initial value of the variable
opts.Upper = [10 15 -3];
opts.StartPoint = [5 5 -9];


%% 
[odorIdxList, ornIdxList] = find(maskPlot == 1);    %find the index of good data

for i = 1 : length(odorIdxList)
    odIdx = odorIdxList(i);     colNum = ornIdxList(i);
    
    rowList = find(input.indOdor == odIdx);
    concPool = input.concList(rowList);
    rPool = input.dff(rowList);
    expIDPool = input.expID(rowList);
    [C, ~, ic] = unique(expIDPool);
    for j = 1 : length(C)
        expList = find(ic == j);
        conc = concPool(expList);
        r = rPool(expList);
        
        
        
        
        
    end
    
    
    
end
















