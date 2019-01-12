function [colorMapSort] = FindColorMap()

load(fullfile('..', 'Figure4','data', 'log10EC50.mat'));

log10EC50(isnan(log10EC50)) = 0;
[~, score, ~,~,~, ~] = pca(log10EC50);
[~, I] = sort(score(:, 1));

myColorMap = jet(34);
colorMapSort = zeros(34, 3);
for i =1:34
    colorMapSort(I(i), :) = myColorMap(i, :);
end

end