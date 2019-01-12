%% load data
load(fullfile('.', 'data', 'dataFigureS3B.mat'));

%% plot
dffSort = dff(idx,:);
ORNListSort = ORNList(idx);

dffRange(:, 1) = min(dffSort, [], 2);
dffRange(:, 2) = max(dffSort, [], 2);

dffStep = dffRange(:, 2) - dffRange(:, 1);
dffGrid = cumsum(dffStep(end: -1: 1));
dffGrid = [0; dffGrid];

figure; set(gcf, 'Position', [418 9 1139 988]);
for i = 21:-1:1
    dffCurrent = dffSort(i, :)+dffGrid(22-i);
    
    plot(t, dffCurrent, 'k'); hold on;
end

axis tight;
yticks(dffGrid(1:21)); yticklabels(ORNListSort(end:-1:1));
xlabel('Time(s)'); 

for i = 1:length(odorT)
    y = [dffGrid(1), dffGrid(end)];
    x = [odorT(i), odorT(i)];
    plot(x, y, '--k');
    
    if mod(i, 2) == 0 && i ~=44
        t = text(odorT(i), dffGrid(end), odorList(i/2));
        t.Rotation = 90;
    end
    
end

% save figure
savefig(gcf, fullfile('results', 'figures', 'traces_FigureS3B.fig'));