function typeMask = SetCurveType(rowNum, colNum)
% Set the type of each dose-response curve
% 1: saturated curves
% 2: highest DF/F value >> 1, and the 2nd highest concentration reach the half maximum
% 3: weak responses, not 1 nor 2, and non-zero.

typeMask = zeros(rowNum, colNum);

typeMask(1, 1) = 3;     %1-pentanol 33b; 
typeMask(1, 2) = 3;     %1-pentanol 45a; 
typeMask(1, 4) = 1;     %1-pentanol 35a; S
typeMask(1, 11) = 3;    %1-pentanol 67b; 
typeMask(1, 13) = 3;    %1-pentanol 13a; 

typeMask(2, 1) = 3;     %3-pentanol 33b;
typeMask(2, 2) = 3;     %3-pentanol 45a;
typeMask(2, 5) = 2;     %3-pentanol 42a;
typeMask(2, 6) = 3;     %3-pentanol 59a;
typeMask(2, 8) = 3;     %3-pentanol 45b;
typeMask(2, 10) = 3;    %3-pentanol 24a;
typeMask(2, 11) = 3;    %3-pentanol 67b; one trial responds
typeMask(2, 12) = 3;    %3-pentanol 85c;
typeMask(2, 13) = 3;    %3-pentanol 13a;
typeMask(2, 16) = 3;    %3-pentanol 22c;
typeMask(2, 17) = 3;    %3-pentanol 42b; one trial responds
typeMask(2, 21) = 3;    %3-pentanol 94a;

typeMask(3, 2) = 3;     %6-methyl-5-hepten-2-ol 45a;
typeMask(3, 4) = 2;     %6-methyl-5-hepten-2-ol 35a;
typeMask(3, 7) = 3;     %6-methyl-5-hepten-2-ol 1a;
typeMask(3, 12) = 1;    %6-methyl-5-hepten-2-ol 85c; S
typeMask(3, 13) = 1;    %6-methyl-5-hepten-2-ol 13a; S

typeMask(4, 1)  = 2;    %3-octanol 33b; 
typeMask(4, 2)  = 3;    %3-octanol 45a; 
typeMask(4, 4)  = 3;    %3-octanol 35a; 
typeMask(4, 7)  = 3;    %3-octanol 1a; 
typeMask(4, 10)  = 3;   %3-octanol 24a; 
typeMask(4, 11)  = 3;   %3-octanol 67b; 
typeMask(4, 12) = 1;    %3-octanol 85c; S
typeMask(4, 13) = 1;    %3-octanol 13a; S

typeMask(5, 1) = 3;     %trans-3-hexen-1-ol 33b; 
typeMask(5, 2) = 3;     %trans-3-hexen-1-ol 45a; 
typeMask(5, 4) = 1;     %trans-3-hexen-1-ol 35a; S
typeMask(5, 5) = 3;     %trans-3-hexen-1-ol 42a; 
typeMask(5, 11) = 2;    %trans-3-hexen-1-ol 67b; 
typeMask(5, 12) = 3;    %trans-3-hexen-1-ol 85c; one trial responds
typeMask(5, 13) = 2;    %trans-3-hexen-1-ol 13a;

typeMask(6, 1) = 3;     %methyl phenyl sulfide 33b;  
typeMask(6, 6) = 3;     %methyl phenyl sulfide 59a;
typeMask(6, 7) = 3;     %methyl phenyl sulfide 1a;
typeMask(6, 8) = 1;     %methyl phenyl sulfide 45b; S
typeMask(6, 10) = 1;     %methyl phenyl sulfide 24a; S
typeMask(6, 11) = 3;    %methyl phenyl sulfide 67b;
typeMask(6, 14) = 2;    %methyl phenyl sulfide 30a;
typeMask(6, 16) = 1;    %methyl phenyl sulfide 22c; S
typeMask(6, 21) = 3;    %methyl phenyl sulfide 94a; 

typeMask(7, 1)  = 3;    %anisole 33b; 
typeMask(7, 2)  = 3;    %anisole 45a; 
typeMask(7, 4)  = 3;    %anisole 35a; one trial responds
typeMask(7, 6)  = 3;    %anisole 59a; 
typeMask(7, 7)  = 3;    %anisole 1a; 
typeMask(7, 8)  = 2;    %anisole 45b; 
typeMask(7, 10) = 1;    %anisole 24a; S
typeMask(7, 11) = 3;    %anisole 67b; 
typeMask(7, 14) = 2;    %anisole 30a; 
typeMask(7, 16) = 2;    %anisole 22c; 
typeMask(7, 21) = 3;    %anisole 94a; 

typeMask(8, 1) = 3;     %2-acetylpyridine 33b;
typeMask(8, 4) = 3;     %2-acetylpyridine 35a;
typeMask(8, 6) = 2;     %2-acetylpyridine 59a;
typeMask(8, 8) = 2;     %2-acetylpyridine 45b;
typeMask(8, 10) = 1;    %2-acetylpyridine 24a; S
typeMask(8, 16) = 2;    %2-acetylpyridine 22c;

typeMask(9, 1) = 3;     %2,5-dimethylpyrazine 33b;
typeMask(9, 2) = 3;     %2,5-dimethylpyrazine 45a; one trial responds
typeMask(9, 6) = 2;     %2,5-dimethylpyrazine 59a;
typeMask(9, 8) = 3;     %2,5-dimethylpyrazine 45b;
typeMask(9, 10)= 3;     %2,5-dimethylpyrazine 67b; one trial responds
typeMask(9, 14) = 3;    %2,5-dimethylpyrazine 30a;
typeMask(9, 16) = 3;    %2,5-dimethylpyrazine 22c;
typeMask(9, 17) = 3;    %2,5-dimethylpyrazine 42b;
typeMask(9, 21) = 3;    %2,5-dimethylpyrazine 94a;

typeMask(10, 1) = 1;    %pentyl acetate 33b; S
typeMask(10, 2) = 2;    %pentyl acetate 45a; not flat yet
typeMask(10, 4) = 1;    %pentyl acetate 35a; S
typeMask(10, 5) = 3;    %pentyl acetate 42a; one trial responds 
typeMask(10, 7) = 2;    %pentyl acetate 1a; not flat yet
typeMask(10, 10) = 2;   %pentyl acetate 24a; 
typeMask(10, 11) = 2;   %pentyl acetate 67b; not flat yet
typeMask(10, 12) = 1;   %pentyl acetate 85c; S
typeMask(10, 13) = 2;   %pentyl acetate 13a;
typeMask(10, 14) = 3;   %pentyl acetate 30a;
typeMask(10, 15) = 3;   %pentyl acetate 82a;
typeMask(10, 16) = 3;   %pentyl acetate 22c;
typeMask(10, 17) = 3;   %pentyl acetate 42b; <0.5

typeMask(11, 1) = 3;    %geranyl acetate 33b;
typeMask(11, 2) = 1;    %geranyl acetate 45a; S
typeMask(11, 7) = 3;    %geranyl acetate 1a; one trial responds
typeMask(11, 10) = 3;    %geranyl acetate 24a;
typeMask(11, 15) = 1;   %geranyl acetate 82a; S
typeMask(11, 17) = 3;   %geranyl acetate 42b;

typeMask(12, 1)  = 3;   %2-methoxyphenyl acetate 33b; 
typeMask(12, 6)  = 2;   %2-methoxyphenyl acetate 59a; 
typeMask(12, 8)  = 3;   %2-methoxyphenyl acetate 45b; 
typeMask(12, 10)  = 3;   %2-methoxyphenyl acetate 24a; 
typeMask(12, 12) = 3;   %2-methoxyphenyl acetate 85c; one trial responds
typeMask(12, 16) = 3;   %2-methoxyphenyl acetate 22c; one trial responds
typeMask(12, 21) = 2;   %2-methoxyphenyl acetate 94a; 

typeMask(13, 1) = 3;    %trans,trans-2,4-nonadienal 33b-47a;  two trial respond
typeMask(13, 2) = 2;    %trans,trans-2,4-nonadienal 45a; 
typeMask(13, 4) = 2;    %trans,trans-2,4-nonadienal 35a; 
typeMask(13, 12) = 3;   %trans,trans-2,4-nonadienal 85c; <0.5 
typeMask(13, 13) = 3;   %trans,trans-2,4-nonadienal 13a; 
typeMask(13, 20) = 1;   %trans,trans-2,4-nonadienal 74a; S

typeMask(14, 3) = 2;    %4m5v 83a;
typeMask(14, 6) = 1;    %4m5v 59a; S
typeMask(14, 8) = 1;    %4m5v 45b; S
typeMask(14, 10) = 3;   %4m5v 24a;
typeMask(14, 14) = 3;   %4m5v 30a;
typeMask(14, 16) = 3;   %4m5v 22c;
typeMask(14, 21) = 3;   %4m5v 94a;

typeMask(15, 3) = 3;    %4,5-dimethylthiazole 83a;
typeMask(15, 6) = 2;    %4,5-dimethylthiazole 59a; 
typeMask(15, 8) = 3;    %4,5-dimethylthiazole 45b;
typeMask(15, 11) = 3;   %4,5-dimethylthiazole 67b;

typeMask(16, 1) = 2;    %4-hexen-3-one 33b;
typeMask(16, 2) = 3;    %4-hexen-3-one 45a;
typeMask(16, 4) = 3;    %4-hexen-3-one 35a;
typeMask(16, 5) = 1;    %4-hexen-3-one 42a; S
typeMask(16, 10) = 3;    %4-hexen-3-one 24a; one trial responds
typeMask(16, 12) = 2;   %4-hexen-3-one 85c;
typeMask(16, 17) = 2;   %4-hexen-3-one 42b;
typeMask(16, 20) = 3;   %4-hexen-3-one 74a;

typeMask(17, 1) = 2;    %2-nonanone 33b; 
typeMask(17, 2) = 3;    %2-nonanone 45a; 
typeMask(17, 4) = 2;    %2-nonanone 35a; 
typeMask(17, 8) = 3;    %2-nonanone 45b; 
typeMask(17, 10) = 2;   %2-nonanone 24a; 
typeMask(17, 11) = 3;   %2-nonanone 67b; 
typeMask(17, 12) = 2;   %2-nonanone 85c; 
typeMask(17, 13) = 2;   %2-nonanone 13a; 
typeMask(17, 15) = 3;   %2-nonanone 82a; one trial responds
typeMask(17, 17) = 3;   %2-nonanone 42b; one trial responds
typeMask(17, 20) = 3;   %2-nonanone 74a; 

typeMask(18, 2) = 3;    %acetal 45a;
typeMask(18, 17) = 2;   %acetal 42b; not flat
typeMask(18, 20) = 3;   %acetal 74a;

typeMask(19, 1) = 3;    %2-phenyl ethanol 33b/47a;
typeMask(19, 4) = 3;    %2-phenyl ethanol 35a;
typeMask(19, 10) = 3;    %2-phenyl ethanol 24a;
typeMask(19, 11) = 2;   %2-phenyl ethanol 67b;
typeMask(19, 12) = 3;   %2-phenyl ethanol 85c;

typeMask(20, 1) = 2;	%butyl acetate 33b/47a;
typeMask(20, 2) = 3;	%butyl acetate 45a;
typeMask(20, 3) = 3;	%butyl acetate 83a; one trial responds
typeMask(20, 4) = 3;	%butyl acetate 35a;
typeMask(20, 5) = 3;	%butyl acetate 42a;
typeMask(20, 10) = 3;	%butyl acetate 24a; one trial responds
typeMask(20, 11) = 3;	%butyl acetate 67b;
typeMask(20, 12) = 1;	%butyl acetate 85c; S
typeMask(20, 13) = 2;	%butyl acetate 13a;
typeMask(20, 15) = 3;	%butyl acetate 82a;
typeMask(20, 16) = 3;	%butyl acetate 22c;
typeMask(20, 17) = 3;	%butyl acetate 42b;
typeMask(20, 21) = 3;	%butyl acetate 94a/94b; one trial responds

typeMask(21, 1) = 3;	%ethyl acetate 33b/47a;
typeMask(21, 3) = 3;	%ethyl acetate 83a; one trial responds
typeMask(21, 4) = 3;	%ethyl acetate 35a;
typeMask(21, 5) = 2;	%ethyl acetate 42a;
typeMask(21, 12) = 3;	%ethyl acetate 85c; one trial responds
typeMask(21, 17) = 1;	%ethyl acetate 42b; S

typeMask(22, 1) = 3;	%benzaldehyde 33b/47a; one trial responds
typeMask(22, 3) = 3;	%benzaldehyde 83a;
typeMask(22, 4) = 2;	%benzaldehyde 35a;
typeMask(22, 5) = 3;	%benzaldehyde 42a; one trial responds
typeMask(22, 7) = 3;	%benzaldehyde 1a;
typeMask(22, 8) = 1;	%benzaldehyde 45b; S
typeMask(22, 9) = 2;	%benzaldehyde 63a;
typeMask(22, 10) = 1;	%benzaldehyde 24a; S
typeMask(22, 11) = 3;	%benzaldehyde 67b;
typeMask(22, 13) = 3;	%benzaldehyde 13a; one trial responds
typeMask(22, 14) = 2;	%benzaldehyde 30a;
typeMask(22, 15) = 3;	%benzaldehyde 82a; one trial responds
typeMask(22, 16) = 3;	%benzaldehydee 22c;
typeMask(22, 20) = 3;	%benzaldehyde 74a;
typeMask(22, 21) = 3;	%benzaldehyde 94a/94b; one trial responds

typeMask(23, 1) = 1;	%2-heptanone 33b/47a; S
typeMask(23, 2) = 2;	%2-heptanone 45a; 
typeMask(23, 3) = 3;	%2-heptanone 83a;  one trial responds
typeMask(23, 4) = 1;	%2-heptanone 35a;  S
typeMask(23, 5) = 1;	%2-heptanone 42a;  S
typeMask(23, 7) = 3;	%2-heptanone 1a;  one trial responds
typeMask(23, 10) = 3;	%2-heptanone 24a; one trial responds
typeMask(23, 11) = 2;	%2-heptanone 67b; 
typeMask(23, 12) = 1;	%2-heptanone 85c; 
typeMask(23, 13) = 2;	%2-heptanone 13a; 
typeMask(23, 15) = 3;	%2-heptanone 82a; one trial responds
typeMask(23, 16) = 3;	%2-heptanone 22c; one trial responds
typeMask(23, 20) = 3;	%2-heptanone 74a; 

typeMask(24, 1) = 3;	%methyl salicylate 33b/47a;
typeMask(24, 6) = 3;	%methyl salicylate 59a;
typeMask(24, 7) = 2;	%methyl salicylate 1a;
typeMask(24, 8) = 3;	%methyl salicylate 45b; one trial responds
typeMask(24, 9) = 2;	%methyl salicylate 63a;
typeMask(24, 10) = 1;	%methyl salicylate 24a; S
typeMask(24, 12) = 2;	%methyl salicylate 85c;
typeMask(24, 13) = 3;	%methyl salicylate 13a; one trial responds
typeMask(24, 16) = 1;	%methyl salicylate 22c; S
typeMask(24, 21) = 3;	%methyl salicylate 94a/94b;

typeMask(25, 1) = 2;	%ethyl butyrate 33b/47a;
typeMask(25, 2) = 3;	%ethyl butyrate 45a;
typeMask(25, 3) = 3;	%ethyl butyrate 83a;
typeMask(25, 4) = 2;	%ethyl butyrate 35a;
typeMask(25, 5) = 2;	%ethyl butyrate 42a;
typeMask(25, 6) = 3;	%ethyl butyrate 59a; one trial responds
typeMask(25, 7) = 3;	%ethyl butyrate 1a; one trial responds
typeMask(25, 9) = 3;	%ethyl butyrate 63a; one trial responds
typeMask(25, 10) = 3;	%ethyl butyrate 24a; one trial responds
typeMask(25, 12) = 3;	%ethyl butyrate 85c;
typeMask(25, 14) = 3;	%ethyl butyrate 30a; one trial responds
typeMask(25, 15) = 3;	%ethyl butyrate 82a;
typeMask(25, 16) = 3;	%ethyl butyrate 22c;
typeMask(25, 17) = 2;	%ethyl butyrate 42b;
typeMask(25, 18) = 2;	%ethyl butyrate 33a;

typeMask(26, 1) = 1;	%isoamyl acetate 33b/47a; S
typeMask(26, 2) = 2;	%isoamyl acetate 45a;
typeMask(26, 3) = 3;	%isoamyl acetate 83a; one trial responds
typeMask(26, 4) = 3;	%isoamyl acetate 35a;
typeMask(26, 5) = 3;	%isoamyl acetate 42a;
typeMask(26, 10) = 3;	%isoamyl acetate 24a; one trial responds
typeMask(26, 11) = 3;	%isoamyl acetate 67b;
typeMask(26, 12) = 1;	%isoamyl acetate 85c; S
typeMask(26, 13) = 3;	%isoamyl acetate 13a;
typeMask(26, 16) = 3;	%isoamyl acetate 22c;
typeMask(26, 17) = 3;	%isoamyl acetate 42b; one trial responds
typeMask(26, 20) = 3;	%isoamyl acetate 74a;

typeMask(27, 1) = 3;	%4-methylcyclohexane 33b/47a; one trial responds
typeMask(27, 3) = 3;	%4-methylcyclohexane 83a;
typeMask(27, 4) = 2;	%4-methylcyclohexane 35a;
typeMask(27, 11) = 2;	%4-methylcyclohexane 67b;
typeMask(27, 12) = 3;	%4-methylcyclohexane 85c; one trial responds
typeMask(27, 16) = 3;	%4-methylcyclohexane 22c;
typeMask(27, 20) = 3;	%4-methylcyclohexane 74a; one trial responds

typeMask(28, 1) = 1;	%hexyl acetate 33b/47a;  S
typeMask(28, 2) = 1;	%hexyl acetate 45a;  S
typeMask(28, 4) = 1;	%hexyl acetate 35a;  S
typeMask(28, 5) = 3;	%hexyl acetate 42a;  
typeMask(28, 7) = 3;	%hexyl acetate 1a;  
typeMask(28, 9) = 3;	%hexyl acetate 63a; one trial responds 
typeMask(28, 11) = 2;	%hexyl acetate 67b;  
typeMask(28, 12) = 3;	%hexyl acetate 85c;  
typeMask(28, 13) = 2;	%hexyl acetate 13a;  
typeMask(28, 15) = 3;	%hexyl acetate 82a;  one trial responds 
typeMask(28, 16) = 3;	%hexyl acetate 22c; 
typeMask(28, 17) = 3;	%hexyl acetate 42b; one trial responds 

typeMask(29, 2) = 3;	%linalool 45a; one trial responds 
typeMask(29, 3) = 3;	%linalool 83a; one trial responds 
typeMask(29, 4) = 3;	%linalool 35a; one trial responds 
typeMask(29, 5) = 3;	%linalool 42a; 
typeMask(29, 6) = 3;	%linalool 59a; 
typeMask(29, 7) = 3;	%linalool 1a; 
typeMask(29, 8) = 3;	%linalool 45b; one trial responds 
typeMask(29, 9) = 3;	%linalool 63a; one trial responds 
typeMask(29, 10) = 3;	%linalool 24a; one trial responds 
typeMask(29, 12) = 2;	%linalool 85c 
typeMask(29, 13) = 2;	%linalool 13a 
typeMask(29, 15) = 3;	%linalool 82a one trial responds 
typeMask(29, 16) = 3;	%linalool 22c one trial responds 
typeMask(29, 17) = 3;	%linalool 42b 
typeMask(29, 21) = 3;	%linalool 94a/94b 

typeMask(30, 1) = 2;	%benzyl acetate 33b/47a; 
typeMask(30, 2) = 2;	%benzyl acetate 45a; 
typeMask(30, 5) = 3;	%benzyl acetate 42a; 
typeMask(30, 6) = 3;	%benzyl acetate 59a; one trial responds 
typeMask(30, 7) = 3;	%benzyl acetate 1a; 
typeMask(30, 9) = 3;	%benzyl acetate 63a; 
typeMask(30,11) = 3;	%benzyl acetate 67b; 
typeMask(30,12) = 3;	%benzyl acetate 85c; 
typeMask(30, 13) = 3;	%benzyl acetate 13a; 
typeMask(30, 15) = 3;	%benzyl acetate 82a; 
typeMask(30, 16) = 3;	%benzyl acetate 22c; 

typeMask(31, 2) = 3;	%4-pheny-2-butanol 45a; one trial responds 
typeMask(31, 4) = 3;	%4-pheny-2-butanol 35a; 
typeMask(31, 7) = 3;	%4-pheny-2-butanol 1a; 
typeMask(31, 8) = 3;	%4-pheny-2-butanol 45b; 
typeMask(31, 9) = 2;	%4-pheny-2-butanol 63a; 
typeMask(31, 10) = 2;	%4-pheny-2-butanol 24a; 
typeMask(31, 11) = 3;	%4-pheny-2-butanol 67b; 
typeMask(31, 12) = 3;	%4-pheny-2-butanol 85c; 
typeMask(31, 13) = 2;	%4-pheny-2-butanol 13a; 
typeMask(31, 14) = 3;	%4-pheny-2-butanol 30a; one trial responds 
typeMask(31, 16) = 3;	%4-pheny-2-butanol 22c; one trial responds 

typeMask(32, 3) = 3;	%myrtenal 83a; 
typeMask(32, 4) = 2;	%myrtenal 35a; 
typeMask(32, 5) = 3;	%myrtenal 42a; one trial responds 
typeMask(32, 8) = 3;	%myrtenal 45b; one trial responds 
typeMask(32, 9) = 2;	%myrtenal 63a; 
typeMask(32, 10) = 2;	%myrtenal 24a; 
typeMask(32, 14) = 3;	%myrtenal 30a; one trial responds 
typeMask(32, 17) = 3;	%myrtenal 42b; one trial responds 
typeMask(32, 19) = 3;	%myrtenal 49a; one trial responds 

typeMask(33, 1) = 3;	%menthol 33b/47a; one trial responds  
typeMask(33, 6) = 3;	%menthol 59a; one trial responds 
typeMask(33, 7) = 2;	%menthol 1a; 
typeMask(33, 8) = 3;	%menthol 45b; 
typeMask(33, 9) = 2;	%menthol 63a; 
typeMask(33, 10) = 3;	%menthol 24a; 
typeMask(33, 16) = 3;	%menthol 22c; one trial responds 
typeMask(33, 17) = 3;	%menthol 42b; 
typeMask(33, 19) = 2;	%menthol 49a; 
typeMask(33, 20) = 3;	%menthol 74a;
typeMask(33, 21) = 3;	%menthol 94a/94b;

typeMask(34, 4) = 1;	%nonane 35a; S
typeMask(34, 5) = 3;	%nonane 42a;
typeMask(34, 11) = 3;   %nonane 67b;
typeMask(34, 12) = 3;   %nonane 85c;
typeMask(34, 13) = 3;   %nonane 13a;
typeMask(34, 16) = 3;   %nonane 22c; one trial responds 
typeMask(34, 17) = 3;   %nonane 42b; one trial responds 
end