function [ connect, synType, inputs, nList ] = GetCon( excelFileName )
%% read the adjacency matrix of the mPN circuit
% output: connect, the matrix, rows are pre-synapse, cols are post-synapse
% synType, synapse type, ORNs and mPNs are ACh, type 1, pickyLNs are type 2.
% inputs, the index of the input nodes, ORNs
% nList, list of neuron IDs.

%read the matrix
if nargin == 0
    excelFileName = fullfile('..', 'data', 'ORN-PickyLN-mPN_Left_Merged.xlsx');
end
connect = xlsread(excelFileName);

ORN = 1:21;
pickyLN = 22:26;
mPN = 27:length(connect(:,1));

inputs = ORN;

synType = zeros(length(connect(:,1)), 1);
synType(ORN) = 1;
synType(pickyLN) = 2;
synType(mPN) = 1;

ORNStr = {'1a ORN left';
'13a ORN left';
'22c ORN left';
'24a ORN left';
'30a ORN left';
'33a ORN left';
'35a ORN left';
'42a ORN left';
'42b ORN left';
'45a ORN left';
'45b ORN left';
'47a & 33b ORN left';
'49a ORN left';
'59a ORN left';
'63a ORN left';
'67b ORN left';
'74a ORN left';
'82a ORN left';
'83a ORN left';
'85c ORN left';
'94a & 94b ORN left';};

pLNStr = {'Picky 0 left';
'Picky 1 left';
'Picky 2 left';
'Picky 3 left';
'Picky 4 left';};

mPNStr = {'mPN iACT A1 left';
'mPN iACT A2 left';
'mPN iACT A3 left';
'mPN iACT B1 left';
'mPN iACT B2 left';
'mPN iACT B3 left';
'mPN iACT C1 left';
'mPN iACT C2 left';
'mPN iACT bilateral LOWER';
'mPN iACT bilateral UPPER';
'mPN iACT VUM';
'mPN D';
'mPN cobra left';
'mPN seahorse left';};

nList = [ORNStr; pLNStr; mPNStr];

end