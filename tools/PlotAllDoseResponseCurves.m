load(fullfile('.', 'data', 'AveRawDataMatrix2ndRound.mat'));

%%
% for i = 20:length(odorList)
for i = 1:length(odorList)
    for j = 1:length(ORNList)
        resp = dataMean(i, j, :);
        if sum(resp) ~= 0
            resp_sem = dataSEM(i, j, :);
            figure;
            errorbar(concList, resp(:), resp_sem(:)); set(gca,'XScale','log');
            title([odorList{i}, '-',num2str(i), '-', ORNList{j}, '-', num2str(j) ]);
            pause(0.01);
        end
    end
end