rawDataTable = dataT;

odorList = input.Odor;
ORNList = input.ORN;

conc = concTs;
dffMean = rspTs;
dffSEM = semTs;

save('ORN_data.mat', 'rawDataTable', 'odorList', 'ORNList', 'conc', ...
    'dffMean', 'dffSEM');