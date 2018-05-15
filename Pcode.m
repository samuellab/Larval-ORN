load('AnalysisResults/FitDoseResponse.mat', 'cMatrix');
load('data/AveRawDataMatrix.mat', 'odorList');
kdM = cMatrix; clear cMatrix

for i = 1:5 %size of primacy code
    % generate primacy code
    pCode = zeros(19, i);
    for j = 1:19
        a = kdM(j, :);
        [b, index] = sort(a, 'ascend');
        temp = index(1:i);
        temp(find(isnan(b(1:i)))) = NaN;
        pCode(j, :) = temp;
    end
    
    % compare cfode
    confuM = zeros(19,19);
    for j = 1:19
        for k = 1:19
            a = pCode(j, :);
            b = pCode(k, :);
            % if a and b have the same elements
            confuM(j, k) = isequaln(sort(a), sort(b));
        end
    end
    figure; imagesc(confuM); title(['Odor identical test using primacy code, p=', num2str(i)]);
    yticks(1:1:19); yticklabels(odorList);
    xticks(1:1:19); xticklabels(odorList);xtickangle(90);
end