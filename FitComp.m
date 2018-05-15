clear; clc;
diary(fullfile('.', 'AnalysisResults', 'DoseResponseFit_log.txt')); diary on;
%% load data, averaged acorss trials of the raw data
load(fullfile('.', 'data', 'AveRawDataMatrix.mat'));
[rowTotal, colTotal, ~] = size(dataRawAve);


%% Label the curves
% 1: saturated curves
% 2: highest DF/F value >> 1, and the 2nd highest concentration reach the half maximum
% 3: weak responses, not 1 nor 2, and non-zero.

maskPlot = zeros(rowTotal, colTotal);

maskPlot(1, 1) = 3;     %1-pentanol 33b; 
maskPlot(1, 2) = 3;     %1-pentanol 45a; 
maskPlot(1, 4) = 1;     %1-pentanol 35a; S
maskPlot(1, 10) = 3;    %1-pentanol 67b; 
maskPlot(1, 11) = 3;    %1-pentanol 85c; 
maskPlot(1, 12) = 3;    %1-pentanol 13a; 

maskPlot(2, 1) = 3;     %3-pentanol 33b;
maskPlot(2, 2) = 3;     %3-pentanol 45a;
maskPlot(2, 5) = 2;     %3-pentanol 42a;
maskPlot(2, 6) = 3;     %3-pentanol 59a;
maskPlot(2, 8) = 3;     %3-pentanol 45b;
maskPlot(2, 9) = 3;     %3-pentanol 24a;
maskPlot(2, 10) = 3;    %3-pentanol 67b;
maskPlot(2, 11) = 3;    %3-pentanol 85c;
maskPlot(2, 12) = 3;    %3-pentanol 13a;
maskPlot(2, 15) = 3;    %3-pentanol 22c;
maskPlot(2, 16) = 3;    %3-pentanol 42b;
maskPlot(2, 18) = 3;    %3-pentanol 94a;

maskPlot(3, 2) = 3;     %6-methyl-5-hepten-2-ol 45a;
maskPlot(3, 4) = 3;     %6-methyl-5-hepten-2-ol 35a;
maskPlot(3, 7) = 3;     %6-methyl-5-hepten-2-ol 1a;
maskPlot(3, 11) = 1;    %6-methyl-5-hepten-2-ol 85c; S
maskPlot(3, 12) = 1;    %6-methyl-5-hepten-2-ol 13a; S

maskPlot(4, 1)  = 2;    %3-octanol 33b; 
maskPlot(4, 2)  = 3;    %3-octanol 45a; 
maskPlot(4, 4)  = 3;    %3-octanol 35a; 
maskPlot(4, 7)  = 3;    %3-octanol 1a; 
maskPlot(4, 9)  = 3;    %3-octanol 24a; 
maskPlot(4, 10)  = 3;   %3-octanol 67b; 
maskPlot(4, 11) = 1;    %3-octanol 85c; S
maskPlot(4, 12) = 1;    %3-octanol 13a; S

maskPlot(5, 1) = 3;     %trans-3-hexen-1-ol 33b; 
maskPlot(5, 2) = 3;     %trans-3-hexen-1-ol 45a; 
maskPlot(5, 4) = 1;     %trans-3-hexen-1-ol 35a; S
maskPlot(5, 5) = 3;     %trans-3-hexen-1-ol 42a; 
maskPlot(5, 10) = 3;    %trans-3-hexen-1-ol 67b; 
maskPlot(5, 11) = 3;    %trans-3-hexen-1-ol 85c;  <0.5
maskPlot(5, 12) = 3;    %trans-3-hexen-1-ol 13a;

maskPlot(6, 1) = 3;     %methyl phenyl sulfide 33b;
maskPlot(6, 6) = 3;     %methyl phenyl sulfide 59a;
maskPlot(6, 7) = 3;     %methyl phenyl sulfide 1a;
maskPlot(6, 8) = 1;     %methyl phenyl sulfide 45b; S
maskPlot(6, 9) = 1;     %methyl phenyl sulfide 24a; S
maskPlot(6, 10) = 3;    %methyl phenyl sulfide 67b;
maskPlot(6, 13) = 2;    %methyl phenyl sulfide 30a;
maskPlot(6, 15) = 1;    %methyl phenyl sulfide 22c; S
maskPlot(6, 18) = 3;    %methyl phenyl sulfide 94a; 

maskPlot(7, 1)  = 3;    %anisole 33b; 
maskPlot(7, 2)  = 3;    %anisole 45a; 
maskPlot(7, 4)  = 3;    %anisole 35a; 
maskPlot(7, 6)  = 3;    %anisole 59a; 
maskPlot(7, 7)  = 3;    %anisole 1a; 
maskPlot(7, 8)  = 3;    %anisole 45b; 
maskPlot(7, 9)  = 1;    %anisole 24a; S
maskPlot(7, 10) = 3;    %anisole 67b; 
maskPlot(7, 13) = 2;    %anisole 30a; 
maskPlot(7, 15) = 3;    %anisole 22c; 
maskPlot(7, 18) = 3;    %anisole 94a; 

maskPlot(8, 1) = 3;     %2-acetylpyridine 33b;
maskPlot(8, 4) = 3;     %2-acetylpyridine 35a;
maskPlot(8, 6) = 3;     %2-acetylpyridine 59a;
maskPlot(8, 8) = 3;     %2-acetylpyridine 45b;
maskPlot(8, 9) = 1;     %2-acetylpyridine 24a; S
maskPlot(8, 15) = 3;    %2-acetylpyridine 22c;

maskPlot(9, 1) = 3;     %2,5-dimethylpyrazine 33b;
maskPlot(9, 2) = 3;     %2,5-dimethylpyrazine 45a;
maskPlot(9, 6) = 3;     %2,5-dimethylpyrazine 59a;
maskPlot(9, 8) = 3;     %2,5-dimethylpyrazine 45b;
maskPlot(9, 10)= 3;     %2,5-dimethylpyrazine 67b;
maskPlot(9, 13) = 3;    %2,5-dimethylpyrazine 30a;
maskPlot(9, 15) = 3;    %2,5-dimethylpyrazine 22c;
maskPlot(9, 16) = 3;    %2,5-dimethylpyrazine 42b;
maskPlot(9, 18) = 3;    %2,5-dimethylpyrazine 94a;

maskPlot(10, 1) = 1;    %pentyl acetate 33b; S
maskPlot(10, 2) = 2;    %pentyl acetate 45a; not flat yet
maskPlot(10, 4) = 1;    %pentyl acetate 35a; S
maskPlot(10, 5) = 3;    %pentyl acetate 42a; 
maskPlot(10, 7) = 2;    %pentyl acetate 1a; not flat yet
maskPlot(10, 9) = 3;    %pentyl acetate 24a; 
maskPlot(10, 10) = 2;   %pentyl acetate 67b; not flat yet
maskPlot(10, 11) = 1;   %pentyl acetate 85c; S
maskPlot(10, 12) = 2;   %pentyl acetate 13a;
maskPlot(10, 13) = 3;   %pentyl acetate 30a;
maskPlot(10, 14) = 3;   %pentyl acetate 82a;
maskPlot(10, 15) = 3;   %pentyl acetate 22c;
maskPlot(10, 16) = 3;   %pentyl acetate 42b;

maskPlot(11, 1) = 1;    %geranyl acetate 33b; S
maskPlot(11, 2) = 1;    %geranyl acetate 45a; S
maskPlot(11, 7) = 3;    %geranyl acetate 1a;
maskPlot(11, 9) = 3;    %geranyl acetate 24a;
maskPlot(11, 14) = 1;   %geranyl acetate 82a; S
maskPlot(11, 16) = 3;   %geranyl acetate 42b;

maskPlot(12, 1)  = 3;   %2-methoxyphenyl acetate 33b; 
maskPlot(12, 6)  = 3;   %2-methoxyphenyl acetate 59a; 
maskPlot(12, 8)  = 3;   %2-methoxyphenyl acetate 45b; 
maskPlot(12, 9)  = 3;   %2-methoxyphenyl acetate 24a; 
maskPlot(12, 11) = 3;   %2-methoxyphenyl acetate 85c; 
maskPlot(12, 15) = 3;   %2-methoxyphenyl acetate 22c; 
maskPlot(12, 18) = 2;   %2-methoxyphenyl acetate 94a; 

maskPlot(13, 2) = 2;    %trans,trans-2,4-nonadienal 45a; 
maskPlot(13, 4) = 3;    %trans,trans-2,4-nonadienal 35a; 
maskPlot(13, 11) = 3;   %trans,trans-2,4-nonadienal 85c; <0.5 
maskPlot(13, 12) = 3;   %trans,trans-2,4-nonadienal 13a; 
maskPlot(13, 17) = 1;   %trans,trans-2,4-nonadienal 74a; S

maskPlot(14, 3) = 2;    %4m5v 83a;
maskPlot(14, 6) = 1;    %4m5v 59a; S
maskPlot(14, 8) = 1;    %4m5v 45b; S
maskPlot(14, 9) = 3;    %4m5v 24a;
maskPlot(14, 15) = 3;   %4m5v 22c;
maskPlot(14, 18) = 3;   %4m5v 94a;
maskPlot(14, 13) = 3;   %4m5v 30a;

maskPlot(15, 3) = 3;    %4,5-dimethylthiazole 83a;
maskPlot(15, 6) = 2;    %4,5-dimethylthiazole 59a; not flat
maskPlot(15, 8) = 3;    %4,5-dimethylthiazole 45b;
maskPlot(15, 10) = 3;   %4,5-dimethylthiazole 67b;

maskPlot(16, 1) = 3;    %4-hexen-3-one 33b;
maskPlot(16, 2) = 3;    %4-hexen-3-one 45a;
maskPlot(16, 4) = 3;    %4-hexen-3-one 35a;
maskPlot(16, 5) = 1;    %4-hexen-3-one 42a; S
maskPlot(16, 9) = 3;    %4-hexen-3-one 24a;
maskPlot(16, 11) = 3;   %4-hexen-3-one 85c;
maskPlot(16, 16) = 2;   %4-hexen-3-one 42b;
maskPlot(16, 17) = 3;   %4-hexen-3-one 74a;

maskPlot(17, 1) = 2;    %2-nonanone 33b; 
maskPlot(17, 2) = 3;    %2-nonanone 45a; 
maskPlot(17, 4) = 3;    %2-nonanone 35a; 
maskPlot(17, 8) = 3;    %2-nonanone 45b; 
maskPlot(17, 9) = 2;    %2-nonanone 24a; 
maskPlot(17, 10) = 3;   %2-nonanone 67b; 
maskPlot(17, 11) = 2;   %2-nonanone 85c; 
maskPlot(17, 12) = 2;   %2-nonanone 13a; 
maskPlot(17, 16) = 3;   %2-nonanone 42b; 
maskPlot(17, 17) = 3;   %2-nonanone 74a; 

maskPlot(18, 2) = 3;    %acetal 45a;
maskPlot(18, 16) = 2;   %acetal 42b; not flat
maskPlot(18, 17) = 3;   %acetal 74a;

maskPlot(19, 4) = 3;    %2-phenyl ethanol 35a;
maskPlot(19, 9) = 3;    %2-phenyl ethanol 24a;
maskPlot(19, 10) = 2;   %2-phenyl ethanol 67b;
maskPlot(19, 11) = 3;   %2-phenyl ethanol 85c;

%% find the non-zero elements
[rows, cols] = find(maskPlot);

%% normalize data using y_max


%% fit using full model
er_full = zeros(1, length(rows));
fun_full = @(x) 1/(1+exp(-b*(x-c)));

a/(1+ exp(-b*(X-c)))


for i = 1:length(rows)
    
end