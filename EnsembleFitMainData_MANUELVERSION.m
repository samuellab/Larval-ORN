addpath(fullfile('.', 'tools'));

%%
fileName = fullfile('data', 'Supplementary Table 3.csv');

dataT = readtable(fileName);    % load Excel file

[input.Odor, ~, input.indOdor] = unique(dataT.Odor, 'stable');	%ORN name
input.ORN = dataT.Properties.VariableNames(4:end);
input.concList = dataT.Concentration;
input.expID = dataT.Exp_ID;
input.dff = table2array(dataT(:, 4:end));

%% average data cross trials for each odor-ORN pair
concTs = zeros(length(input.Odor), length(input.ORN), 5);
rspTs  = concTs;

for i = 1 : length(input.Odor)
    odI = find(input.indOdor == i); % row index of each odor
    
    concPool = input.concList(odI); % concentration list
    dffPool = input.dff(odI, :);    % response data block
    
    [expList, ia, ic] = unique(input.expID(odI), 'stable');   % find out how many trials each odor
    
    dffPoolTensor = reshape(dffPool', length(input.ORN), [], length(expList));
    concPoolMat = reshape(concPool, [], length(expList)); 
    
    for j = 1 : length(input.ORN)
        dffMat = squeeze(dffPoolTensor(j, :, :));
        
        % each colume is a trial, find trilas are not NaN
        colIdx = find(~isnan(dffMat(1, :)));
        
        % check if the concentration list is the same for these trials
        concBlock =concPoolMat(:, colIdx);
        dffBlock = dffMat(:, colIdx);
        if isequal(concBlock(:, 1:end-1), concBlock(:, 2:end))
            concTs(i, j, :) = concBlock(:, 1);
            rspTs(i, j, :) = mean(dffBlock, 2);
        else
            error('Trials do not share the same set of odor concentration.');
        end
    end
end

%%
[rowTotal, colTotal, ~] = size(concTs);
maskPlot = zeros(rowTotal, colTotal);

maskPlot(1, 1) = 3;     %1-pentanol 33b; 
% maskPlot(1, 2) = 3;     %1-pentanol 45a; 
maskPlot(1, 4) = 1;     %1-pentanol 35a; S
maskPlot(1, 11) = 3;    %1-pentanol 67b; 
maskPlot(1, 13) = 3;    %1-pentanol 13a; 

maskPlot(2, 1) = 3;     %3-pentanol 33b;
maskPlot(2, 2) = 3;     %3-pentanol 45a;
maskPlot(2, 5) = 2;     %3-pentanol 42a;
maskPlot(2, 6) = 3;     %3-pentanol 59a;
maskPlot(2, 8) = 3;     %3-pentanol 45b;
maskPlot(2, 10) = 3;     %3-pentanol 24a;
maskPlot(2, 11) = 3;    %3-pentanol 67b;
maskPlot(2, 12) = 3;    %3-pentanol 85c;
maskPlot(2, 13) = 3;    %3-pentanol 13a;
maskPlot(2, 16) = 3;    %3-pentanol 22c;
maskPlot(2, 17) = 3;    %3-pentanol 42b;
maskPlot(2, 21) = 3;    %3-pentanol 94a;

maskPlot(3, 2) = 3;     %6-methyl-5-hepten-2-ol 45a;
maskPlot(3, 4) = 3;     %6-methyl-5-hepten-2-ol 35a;
maskPlot(3, 7) = 3;     %6-methyl-5-hepten-2-ol 1a;
maskPlot(3, 12) = 1;    %6-methyl-5-hepten-2-ol 85c; S
maskPlot(3, 13) = 1;    %6-methyl-5-hepten-2-ol 13a; S

maskPlot(4, 1)  = 2;    %3-octanol 33b; 
maskPlot(4, 2)  = 3;    %3-octanol 45a; 
maskPlot(4, 4)  = 3;    %3-octanol 35a; 
maskPlot(4, 7)  = 3;    %3-octanol 1a; 
maskPlot(4, 10)  = 3;    %3-octanol 24a; 
maskPlot(4, 11)  = 3;   %3-octanol 67b; 
maskPlot(4, 12) = 1;    %3-octanol 85c; S
maskPlot(4, 13) = 1;    %3-octanol 13a; S

maskPlot(5, 1) = 3;     %trans-3-hexen-1-ol 33b; 
maskPlot(5, 2) = 3;     %trans-3-hexen-1-ol 45a; 
maskPlot(5, 4) = 1;     %trans-3-hexen-1-ol 35a; S
maskPlot(5, 5) = 3;     %trans-3-hexen-1-ol 42a; 
maskPlot(5, 11) = 3;    %trans-3-hexen-1-ol 67b; 
maskPlot(5, 12) = 3;    %trans-3-hexen-1-ol 85c;  <0.5
maskPlot(5, 13) = 3;    %trans-3-hexen-1-ol 13a;

maskPlot(6, 1) = 3;     %methyl phenyl sulfide 33b;  
maskPlot(6, 6) = 3;     %methyl phenyl sulfide 59a;
maskPlot(6, 7) = 3;     %methyl phenyl sulfide 1a;
maskPlot(6, 8) = 1;     %methyl phenyl sulfide 45b; S
maskPlot(6, 10) = 1;     %methyl phenyl sulfide 24a; S
maskPlot(6, 11) = 3;    %methyl phenyl sulfide 67b;
maskPlot(6, 14) = 2;    %methyl phenyl sulfide 30a;
maskPlot(6, 16) = 1;    %methyl phenyl sulfide 22c; S
% maskPlot(6, 21) = 3;    %methyl phenyl sulfide 94a; 

maskPlot(7, 1)  = 3;    %anisole 33b; 
maskPlot(7, 2)  = 3;    %anisole 45a; 
maskPlot(7, 4)  = 3;    %anisole 35a; 
maskPlot(7, 6)  = 3;    %anisole 59a; 
maskPlot(7, 7)  = 3;    %anisole 1a; 
maskPlot(7, 8)  = 3;    %anisole 45b; 
maskPlot(7, 10)  = 1;    %anisole 24a; S
maskPlot(7, 11) = 3;    %anisole 67b; 
maskPlot(7, 14) = 2;    %anisole 30a; 
maskPlot(7, 16) = 3;    %anisole 22c; 
maskPlot(7, 21) = 3;    %anisole 94a; 

maskPlot(8, 1) = 3;     %2-acetylpyridine 33b;
maskPlot(8, 4) = 3;     %2-acetylpyridine 35a;
maskPlot(8, 6) = 3;     %2-acetylpyridine 59a;
maskPlot(8, 8) = 3;     %2-acetylpyridine 45b;
maskPlot(8, 10) = 1;     %2-acetylpyridine 24a; S
maskPlot(8, 16) = 3;    %2-acetylpyridine 22c;

maskPlot(9, 1) = 3;     %2,5-dimethylpyrazine 33b;
maskPlot(9, 2) = 3;     %2,5-dimethylpyrazine 45a;
maskPlot(9, 6) = 3;     %2,5-dimethylpyrazine 59a;
maskPlot(9, 8) = 3;     %2,5-dimethylpyrazine 45b;
% maskPlot(9, 10)= 3;     %2,5-dimethylpyrazine 67b; <0.5
maskPlot(9, 14) = 3;    %2,5-dimethylpyrazine 30a;
maskPlot(9, 16) = 3;    %2,5-dimethylpyrazine 22c;
maskPlot(9, 17) = 3;    %2,5-dimethylpyrazine 42b;
maskPlot(9, 21) = 3;    %2,5-dimethylpyrazine 94a;

maskPlot(10, 1) = 1;    %pentyl acetate 33b; S
maskPlot(10, 2) = 2;    %pentyl acetate 45a; not flat yet
maskPlot(10, 4) = 1;    %pentyl acetate 35a; S
maskPlot(10, 5) = 3;    %pentyl acetate 42a;  <0.5
maskPlot(10, 7) = 2;    %pentyl acetate 1a; not flat yet
maskPlot(10, 10) = 3;   %pentyl acetate 24a; 
maskPlot(10, 11) = 2;   %pentyl acetate 67b; not flat yet
maskPlot(10, 12) = 1;   %pentyl acetate 85c; S
maskPlot(10, 13) = 2;   %pentyl acetate 13a;
maskPlot(10, 14) = 3;   %pentyl acetate 30a;
maskPlot(10, 15) = 3;   %pentyl acetate 82a;
maskPlot(10, 16) = 3;   %pentyl acetate 22c;
maskPlot(10, 17) = 3;   %pentyl acetate 42b; <0.5

maskPlot(11, 1) = 1;    %geranyl acetate 33b; S
maskPlot(11, 2) = 1;    %geranyl acetate 45a; S
maskPlot(11, 7) = 3;    %geranyl acetate 1a; <0.5
maskPlot(11, 10) = 3;    %geranyl acetate 24a;
maskPlot(11, 15) = 1;   %geranyl acetate 82a; S
maskPlot(11, 17) = 3;   %geranyl acetate 42b;

maskPlot(12, 1)  = 3;   %2-methoxyphenyl acetate 33b; 
maskPlot(12, 6)  = 3;   %2-methoxyphenyl acetate 59a; 
maskPlot(12, 8)  = 3;   %2-methoxyphenyl acetate 45b; 
maskPlot(12, 10)  = 3;   %2-methoxyphenyl acetate 24a; 
maskPlot(12, 12) = 3;   %2-methoxyphenyl acetate 85c; 
maskPlot(12, 16) = 3;   %2-methoxyphenyl acetate 22c; 
maskPlot(12, 21) = 2;   %2-methoxyphenyl acetate 94a; 

maskPlot(13, 2) = 2;    %trans,trans-2,4-nonadienal 45a; 
maskPlot(13, 4) = 3;    %trans,trans-2,4-nonadienal 35a; 
maskPlot(13, 12) = 3;   %trans,trans-2,4-nonadienal 85c; <0.5 
maskPlot(13, 13) = 3;   %trans,trans-2,4-nonadienal 13a; 
maskPlot(13, 20) = 1;   %trans,trans-2,4-nonadienal 74a; S

maskPlot(14, 3) = 2;    %4m5v 83a;
maskPlot(14, 6) = 1;    %4m5v 59a; S
maskPlot(14, 8) = 1;    %4m5v 45b; S
maskPlot(14, 10) = 3;    %4m5v 24a;
maskPlot(14, 14) = 3;   %4m5v 30a;
maskPlot(14, 16) = 3;   %4m5v 22c;
maskPlot(14, 21) = 3;   %4m5v 94a;

maskPlot(15, 3) = 3;    %4,5-dimethylthiazole 83a;
maskPlot(15, 6) = 2;    %4,5-dimethylthiazole 59a; not flat
maskPlot(15, 8) = 3;    %4,5-dimethylthiazole 45b;
maskPlot(15, 11) = 3;   %4,5-dimethylthiazole 67b;

maskPlot(16, 1) = 3;    %4-hexen-3-one 33b;
maskPlot(16, 2) = 3;    %4-hexen-3-one 45a;
maskPlot(16, 4) = 3;    %4-hexen-3-one 35a;
maskPlot(16, 5) = 1;    %4-hexen-3-one 42a; S
maskPlot(16, 10) = 3;    %4-hexen-3-one 24a;
maskPlot(16, 12) = 3;   %4-hexen-3-one 85c;
maskPlot(16, 17) = 2;   %4-hexen-3-one 42b;
maskPlot(16, 20) = 3;   %4-hexen-3-one 74a;

maskPlot(17, 1) = 2;    %2-nonanone 33b; 
maskPlot(17, 2) = 3;    %2-nonanone 45a; 
maskPlot(17, 4) = 3;    %2-nonanone 35a; 
maskPlot(17, 8) = 3;    %2-nonanone 45b; 
maskPlot(17, 10) = 2;    %2-nonanone 24a; 
maskPlot(17, 11) = 3;   %2-nonanone 67b; 
maskPlot(17, 12) = 2;   %2-nonanone 85c; 
maskPlot(17, 13) = 2;   %2-nonanone 13a; 
maskPlot(17, 17) = 3;   %2-nonanone 42b; 
maskPlot(17, 20) = 3;   %2-nonanone 74a; 

maskPlot(18, 2) = 3;    %acetal 45a;
maskPlot(18, 17) = 2;   %acetal 42b; not flat
maskPlot(18, 20) = 3;   %acetal 74a;

maskPlot(19, 4) = 3;    %2-phenyl ethanol 35a;
maskPlot(19, 10) = 3;    %2-phenyl ethanol 24a;
maskPlot(19, 11) = 2;   %2-phenyl ethanol 67b;
maskPlot(19, 12) = 3;   %2-phenyl ethanol 85c;

maskPlot(20, 1) = 2;	%butyl acetate 33b/47a;
maskPlot(20, 2) = 3;	%butyl acetate 45a;
maskPlot(20, 3) = 3;	%butyl acetate 83a;
maskPlot(20, 4) = 3;	%butyl acetate 35a;
maskPlot(20, 5) = 3;	%butyl acetate 42a;
% maskPlot(20, 10) = 3;	%butyl acetate 24a;
maskPlot(20, 11) = 3;	%butyl acetate 67b;
maskPlot(20, 12) = 1;	%butyl acetate 85c; S
maskPlot(20, 13) = 3;	%butyl acetate 13a;
maskPlot(20, 15) = 3;	%butyl acetate 82a;
maskPlot(20, 16) = 3;	%butyl acetate 22c;
maskPlot(20, 17) = 3;	%butyl acetate 42b;
% maskPlot(20, 21) = 3;	%butyl acetate 94a/94b;

maskPlot(21, 1) = 3;	%ethyl acetate 33b/47a;
% maskPlot(21, 3) = 3;	%ethyl acetate 83a;
% maskPlot(21, 4) = 3;	%ethyl acetate 35a;
maskPlot(21, 5) = 3;	%ethyl acetate 42a;
% maskPlot(21, 12) = 3;	%ethyl acetate 85c;
maskPlot(21, 17) = 1;	%ethyl acetate 42b; S

% maskPlot(22, 1) = 3;	%benzaldehyde 33b/47a;
maskPlot(22, 3) = 3;	%benzaldehyde 83a;
maskPlot(22, 4) = 3;	%benzaldehyde 35a;
maskPlot(22, 5) = 3;	%benzaldehyde 42a;
maskPlot(22, 7) = 3;	%benzaldehyde 1a;
maskPlot(22, 8) = 1;	%benzaldehyde 45b; S
maskPlot(22, 9) = 3;	%benzaldehyde 63a;
maskPlot(22, 10) = 1;	%benzaldehyde 24a; S
maskPlot(22, 11) = 3;	%benzaldehyde 67b;
% maskPlot(22, 13) = 3;	%benzaldehyde 13a;
maskPlot(22, 14) = 3;	%benzaldehyde 30a;
% maskPlot(22, 15) = 3;	%benzaldehyde 82a;
maskPlot(22, 16) = 3;	%benzaldehydee 22c;
% maskPlot(22, 20) = 3;	%benzaldehyde 74a;
% maskPlot(22, 21) = 3;	%benzaldehyde 94a/94b;

maskPlot(23, 1) = 1;	%2-heptanone 33b/47a; S
maskPlot(23, 2) = 3;	%2-heptanone 45a; 
maskPlot(23, 3) = 3;	%2-heptanone 83a; 
maskPlot(23, 4) = 1;	%2-heptanone 35a;  S
maskPlot(23, 5) = 1;	%2-heptanone 42a;  S
maskPlot(23, 7) = 3;	%2-heptanone 1a; 
maskPlot(23, 10) = 3;	%2-heptanone 24a; 
maskPlot(23, 11) = 2;	%2-heptanone 67b; 
maskPlot(23, 12) = 4;	%2-heptanone 85c; 
maskPlot(23, 13) = 2;	%2-heptanone 13a; 
maskPlot(23, 15) = 3;	%2-heptanone 82a; 
% maskPlot(23, 16) = 3;	%2-heptanone 22c; 
maskPlot(23, 20) = 3;	%2-heptanone 74a; 

% maskPlot(24, 1) = 3;	%methyl salicylate 33b/47a;
maskPlot(24, 6) = 3;	%methyl salicylate 59a;
maskPlot(24, 7) = 2;	%methyl salicylate 1a;
maskPlot(24, 8) = 3;	%methyl salicylate 45b;
maskPlot(24, 9) = 2;	%methyl salicylate 63a;
maskPlot(24, 10) = 1;	%methyl salicylate 24a; S
maskPlot(24, 12) = 2;	%methyl salicylate 85c;
% maskPlot(24, 13) = 3;	%methyl salicylate 13a;
maskPlot(24, 16) = 4;	%methyl salicylate 22c; S
% maskPlot(24, 21) = 3;	%methyl salicylate 94a/94b;

maskPlot(25, 1) = 3;	%ethyl butyrate 33b/47a;
maskPlot(25, 2) = 3;	%ethyl butyrate 45a;
maskPlot(25, 3) = 3;	%ethyl butyrate 83a;
maskPlot(25, 4) = 3;	%ethyl butyrate 35a;
maskPlot(25, 5) = 2;	%ethyl butyrate 42a;
maskPlot(25, 6) = 3;	%ethyl butyrate 59a;
% maskPlot(25, 7) = 3;	%ethyl butyrate 1a;
maskPlot(25, 9) = 3;	%ethyl butyrate 63a;
% maskPlot(25, 10) = 3;	%ethyl butyrate 24a;
maskPlot(25, 12) = 3;	%ethyl butyrate 85c;
maskPlot(25, 14) = 3;	%ethyl butyrate 30a;
maskPlot(25, 15) = 3;	%ethyl butyrate 82a;
maskPlot(25, 16) = 3;	%ethyl butyrate 22c;
maskPlot(25, 17) = 2;	%ethyl butyrate 42b;
maskPlot(25, 18) = 3;	%ethyl butyrate 33a;

maskPlot(26, 1) = 1;	%isoamyl acetate 33b/47a; S
maskPlot(26, 2) = 3;	%isoamyl acetate 45a;
maskPlot(26, 3) = 3;	%isoamyl acetate 83a;
maskPlot(26, 4) = 3;	%isoamyl acetate 35a;
maskPlot(26, 5) = 3;	%isoamyl acetate 42a;
% maskPlot(26, 6) = 3;	%isoamyl acetate 59a;
% maskPlot(26, 10) = 3;	%isoamyl acetate 24a;
maskPlot(26, 11) = 3;	%isoamyl acetate 67b;
maskPlot(26, 12) = 1;	%isoamyl acetate 85c; S
maskPlot(26, 13) = 3;	%isoamyl acetate 13a;
maskPlot(26, 16) = 3;	%isoamyl acetate 22c;
maskPlot(26, 17) = 3;	%isoamyl acetate 42b;
maskPlot(26, 20) = 3;	%isoamyl acetate 74a;

maskPlot(27, 1) = 3;	%4-methylcyclohexane 33b/47a;
maskPlot(27, 3) = 3;	%4-methylcyclohexane 83a;
maskPlot(27, 4) = 3;	%4-methylcyclohexane 35a;
maskPlot(27, 11) = 2;	%4-methylcyclohexane 67b;
maskPlot(27, 12) = 3;	%4-methylcyclohexane 85c;
maskPlot(27, 16) = 3;	%4-methylcyclohexane 22c;
maskPlot(27, 20) = 3;	%4-methylcyclohexane 74a;

maskPlot(28, 1) = 1;	%hexyl acetate 33b/47a;  S
maskPlot(28, 2) = 1;	%hexyl acetate 45a;  S
maskPlot(28, 4) = 1;	%hexyl acetate 35a;  S
maskPlot(28, 5) = 3;	%hexyl acetate 42a;  
maskPlot(28, 7) = 3;	%hexyl acetate 1a;  
maskPlot(28, 9) = 3;	%hexyl acetate 63a;  
maskPlot(28, 11) = 2;	%hexyl acetate 67b;  
maskPlot(28, 12) = 2;	%hexyl acetate 85c;  
maskPlot(28, 13) = 2;	%hexyl acetate 13a;  
% maskPlot(28, 15) = 3;	%hexyl acetate 82a;  
maskPlot(28, 16) = 3;	%hexyl acetate 22c; 
% maskPlot(28, 17) = 3;	%hexyl acetate 42b; 

% maskPlot(29, 2) = 3;	%linalool 45a; 
% maskPlot(29, 3) = 3;	%linalool 83a; 
% maskPlot(29, 4) = 3;	%linalool 35a; 
maskPlot(29, 5) = 3;	%linalool 42a; 
maskPlot(29, 6) = 3;	%linalool 59a; 
maskPlot(29, 7) = 3;	%linalool 1a; 
maskPlot(29, 8) = 3;	%linalool 45b; 
maskPlot(29, 9) = 3;	%linalool 63a; 
% maskPlot(29, 10) = 3;	%linalool 24a; 
maskPlot(29, 12) = 2;	%linalool 85c 
maskPlot(29, 13) = 2;	%linalool 13a 
% maskPlot(29, 15) = 3;	%linalool 82a 
maskPlot(29, 16) = 3;	%linalool 22c 
maskPlot(29, 17) = 3;	%linalool 42b 
% maskPlot(29, 21) = 3;	%linalool 94a/94b 

maskPlot(30, 1) = 2;	%benzyl acetate 33b/47a; 
maskPlot(30, 2) = 2;	%benzyl acetate 45a; 
% maskPlot(30, 5) = 3;	%benzyl acetate 42a; 
% maskPlot(30, 6) = 3;	%benzyl acetate 59a; 
maskPlot(30, 7) = 3;	%benzyl acetate 1a; 
maskPlot(30, 9) = 2;	%benzyl acetate 63a; 
maskPlot(30,11) = 3;	%benzyl acetate 67b; 
maskPlot(30,12) = 3;	%benzyl acetate 85c; 
maskPlot(30, 13) = 3;	%benzyl acetate 13a; 
maskPlot(30, 15) = 3;	%benzyl acetate 82a; 
maskPlot(30, 16) = 3;	%benzyl acetate 22c; 

% maskPlot(31, 2) = 3;	%4-pheny-2-butanol 45a; 
maskPlot(31, 4) = 3;	%4-pheny-2-butanol 35a; 
maskPlot(31, 7) = 3;	%4-pheny-2-butanol 1a; 
maskPlot(31, 8) = 3;	%4-pheny-2-butanol 45b; 
maskPlot(31, 9) = 2;	%4-pheny-2-butanol 63a; 
maskPlot(31, 10) = 2;	%4-pheny-2-butanol 24a; 
maskPlot(31, 11) = 3;	%4-pheny-2-butanol 67b; 
maskPlot(31, 12) = 3;	%4-pheny-2-butanol 85c; 
maskPlot(31, 13) = 3;	%4-pheny-2-butanol 13a; 
% maskPlot(31, 14) = 3;	%4-pheny-2-butanol 30a; 
maskPlot(31, 16) = 3;	%4-pheny-2-butanol 22c; 

maskPlot(32, 3) = 3;	%myrtenal 83a; 
maskPlot(32, 4) = 2;	%myrtenal 35a; 
% maskPlot(32, 5) = 3;	%myrtenal 42a; 
maskPlot(32, 8) = 3;	%myrtenal 45b; 
maskPlot(32, 9) = 3;	%myrtenal 63a; 
maskPlot(32, 10) = 3;	%myrtenal 24a; 
% maskPlot(32, 14) = 3;	%myrtenal 30a; 
% maskPlot(32, 17) = 3;	%myrtenal 42b; 
maskPlot(32, 19) = 3;	%myrtenal 49a; 

% maskPlot(33, 1) = 3;	%menthol 33b/47a; 
% maskPlot(33, 6) = 3;	%menthol 59a; 
maskPlot(33, 7) = 3;	%menthol 1a; 
maskPlot(33, 8) = 3;	%menthol 45b; 
maskPlot(33, 9) = 2;	%menthol 63a; 
% maskPlot(33, 10) = 3;	%menthol 24a; 
% maskPlot(33, 16) = 3;	%menthol 22c; 
maskPlot(33, 17) = 3;	%menthol 42b; 
maskPlot(33, 19) = 2;	%menthol 49a; 
maskPlot(33, 20) = 3;	%menthol 74a;
maskPlot(33, 21) = 3;	%menthol 94a/94b;

maskPlot(34, 4) = 1;	%nonane 35a; S
maskPlot(34, 5) = 3;	%nonane 42a;
maskPlot(34, 11) = 3;   %nonane 67b;
maskPlot(34, 12) = 3;   %nonane 85c;
maskPlot(34, 13) = 3;   %nonane 13a;
% maskPlot(34, 16) = 3;   %nonane 22c; <0.1
% maskPlot(34, 17) = 3;   %nonane 42b; <0.5


%% go though the paris find all the non-active paris. 
% fitMask = 99 * ones(length(input.Odor), length(input.ORN));
% 
% rMatSum = sum(rspTs, 3); 
% 
% fitMask(find(rMatSum == 0)) = 0; 
% countPairZeros = length(find(rMatSum(:) == 0));
% disp(['Count of non activated odor-ORN pairs:', num2str(countPairZeros), ...
%     ', which is ',num2str(100*countPairZeros/numel(rMatSum)), ...
%     '% of the total experiments.']);

%% pre-fit, see which pair could be fitted well
% set up the fitting
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [0 0 -11]; % setup the range and initial value of the variable
opts.Upper = [max(rspTs(:))*1.5 10 0];
opts.StartPoint = [4 4 -7];

cMatrix = NaN(rowTotal, colTotal);     % the thresholds for each odor_ORN pair, 'c' in the equation
aMatrix = NaN(rowTotal, colTotal);     % the saturated amplitude for each odor-ORN pair
r2Matrix= NaN(rowTotal, colTotal);     % the R-square value of the fitting

%%
disp('----------PRE-FIT----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftRow, ftCol] = find(maskPlot == 1);
gfX = [];  gfY = [];  
gfRC = [];  gfCoeff = [];	gfR2 = [];
for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft, opts);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
%     flag = coeff(1)>2 && coeff(3) < -5.8;
%     if rSq > 0.9 && flag
%        fitMask(ftRow(i), ftCol(i)) = 1;
       gfX = [gfX; log10(xx')];  gfY = [gfY; yy'];  gfRC = [gfRC; ftRow(i), ftCol(i)]; 
       gfCoeff = [gfCoeff; coeff];  gfR2 = [gfR2; rSq];
       
       fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            input.Odor{gfRC(end, 1)}, input.ORN{gfRC(end, 2)}, ...
            gfCoeff(end, 1), gfCoeff(end, 2), gfCoeff(end, 3),gfR2(end));
        
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);
        
%         figure; plot(log10(xx), yy, 'ok'); hold on;
%         plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
%         title([input.Odor{gfRC(end, 1)}, input.ORN{gfRC(end, 2)}]);
%     end
end

%% ensemble fit of all these good fitted curves
addpath(fullfile('.', 'tools'));

slop0 = median(gfCoeff(:, 2)); ampVec0 = gfCoeff(:, 1); kdVec0 = gfCoeff(:, 3);

disp('----------FIT CURVE ENSEMBLE:----------');
fprintf('%-5s\t%-5s\t\n', 'Slop', 'R^2');
[slop, ampVec, kdVec, rSquare] = EnsembleMiniSearch(gfX, gfY, hillEq, slop0, ampVec0, kdVec0);

fprintf('%.2f\t%.2f\t\n', slop, rSquare);

% plot the data 
dataXEn = gfX -  repmat(kdVec,  1, length(gfX(1,:)));
dataYEn = gfY ./ repmat(ampVec, 1, length(gfY(1,:)));

% save into results
results.fitCoeffFMS = [ampVec, repmat(slop, length(ampVec), 1), kdVec];

figure; 
plot(dataXEn(:), dataYEn(:), 'ok'); hold on;
xPlot = linspace(min(dataXEn(:)), max(dataXEn(:)), 100);
yPlot = hillEq(1, slop, 0, xPlot);
plot(xPlot, yPlot, 'r'); xlabel('Relative log10(c)'); ylabel('Norm.(\DeltaF/F)');
hold off;
for i = 1:length(ampVec)
    cMatrix(gfRC(i,1), gfRC(i,2)) = kdVec(i);
    aMatrix(gfRC(i,1), gfRC(i,2)) = ampVec(i);
    r2Matrix(gfRC(i,1), gfRC(i,2))= gfR2(i);
end


%% fit with fixed slop
ft2 = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Robust = 'Bisquare';
opts2.Lower = [0 slop -11]; % setup the range and initial value of the variable
opts2.Upper = [max(rspTs(:))*1.5 slop 0];
opts2.StartPoint = [4 slop -7];

disp('----------FIT WITH KNOWN SLOP----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftRow, ftCol] = find(maskPlot == 2);
gfX2 = [];  gfY2 = [];
gfRC2 = [];  gfCoeff2 = [];	gfR22 = [];
for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft2, opts2);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
%     flag = coeff(1)>2 && coeff(3) < -5;
%     if rSq > 0.9 && flag
%        fitMask(ftRow(i), ftCol(i)) = 2;
       gfX2 = [gfX2; log10(xx')];  gfY2 = [gfY2; yy'];  gfRC2 = [gfRC2; ftRow(i), ftCol(i)]; 
       gfCoeff2 = [gfCoeff2; coeff];  gfR22 = [gfR22; rSq];
       
       cMatrix(ftRow(i), ftCol(i)) = coeff(3);
       aMatrix(ftRow(i), ftCol(i)) = coeff(1);
       r2Matrix(ftRow(i), ftCol(i))= rSq;
       
       fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            input.Odor{gfRC2(end, 1)}, input.ORN{gfRC2(end, 2)}, ...
            gfCoeff2(end, 1), gfCoeff2(end, 2), gfCoeff2(end, 3),gfR22(end));
        
%         xP = linspace(min(log10(xx)), max(log10(xx)), 50);
%         yP = hillEq(coeff(1), coeff(2), coeff(3), xP);
        
%         figure; 
%         plot(log10(xx), yy, 'ok'); hold on;
%         plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
%         title([input.Odor{gfRC2(end, 1)}, input.ORN{gfRC2(end, 2)}]);
%     end
end


%% Estimate the maximum values for each odor (for fitting the type 3 curves)
maxAofOdor = zeros(length(input.Odor), 1); % define parameter to store maximum amplitude 
lowerBound = 3; % set a lower bound of the maximum amplitude

disp('----------FIND AMPLITUDE OF EACH ODOR FOR TYPE 3 CURVE FITTING----------');
fprintf('%30s\t%-5s\t\n', 'Odor Name', 'Amplitude');

maxProj = max(rspTs, [], 3); %max of the data along the concentration
 
for i = 1:length(input.Odor)   
    index = ~isnan(aMatrix(i, :)); %find the nun-zero elements of each row (odor), from the saturated dataset
    maxSeqTemp = aMatrix(i, index);
    if sum(maxSeqTemp>lowerBound)>0
        maxSeq = maxSeqTemp(maxSeqTemp>lowerBound); %select the max values higher than the lower bound
        maxEst = mean(maxSeq); % use the mean as the maximum
        if maxEst>lowerBound
            maxAofOdor(i) = maxEst;
        else %if there is no qualified elements, find the maximum in the raw data
            indexBackup = find(maxProj(i, :)>lowerBound);
            seqBackup = maxProj(i, indexBackup);
            maxEstBackup = mean(seqBackup);
            maxAofOdor(i) = maxEstBackup;
        end
    else
        
    end

	fprintf('%30s\t%.2f\t\n', input.Odor{i}, maxAofOdor(i));
end

%%
