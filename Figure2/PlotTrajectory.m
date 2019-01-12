load(fullfile('.', 'data', 'dataFigureS3D.mat'));

colormap1 = [204,236,230; ...
153,216,201; ...
102,194,164; ...
44,162,95; ...
0,109,44]/255;

colormap2 = [253,212,158; ...
253,187,132; ...
252,141,89; ...
227,74,51; ...
179,0,0]/255;

figure; hold on;

for ff = 1:2
    
    switch ff
        case 1
            data = data1; colormap = colormap1;
        case 2
            data = data2; colormap = colormap2;
    end
    
dff = data.dff;
t = data.t;
id = data.id;
oTInfo = data.oTInfo;

ps = 2:2:10;
ptStart = oTInfo(ps);
ptEnd = oTInfo(ps+1) + 15;

%% load ORN list
load(fullfile('.', 'results', 'doseResponseData.mat'), 'ORNList');

lut = zeros(1, length(id));
for i = 1:length(id)
    myStr = ['Or', id{i}];
    if ~isempty(find(strcmp(ORNList, myStr)))
        lut(i) = find(strcmp(ORNList, myStr));
    else
%         lut(i) = [];
%         dff(i, :) = [];
    end
end
if ~isempty(find(~lut))
    lut(find(~lut)) = [];
	dff(find(~lut), :) = [];
end

%% 
for i = 1:5
    tStart = ptStart(i); tEnd = ptEnd(i);
    [~, iS]= min(abs(t-tStart)); [~, iE]= min(abs(t-tEnd));
    temp= dff(:, iS-1 : iE);
    
    box = zeros(21, size(temp, 2));
    for j = 1:length(lut)
        box(lut(j), :) = temp(j, :);
    end
    
    dffSeq{i} = box;
end

for i = 1:5
    a = dffSeq{i};
    a = a';
    
    myScore = a/(coeff');
    
    scatter3(myScore(:, 1), myScore(:, 2), myScore(:, 3), ...
        'MarkerFaceColor', colormap(i, :), ...
        'MarkerEdgeColor', 'w'); 
    plot3(myScore(:, 1), myScore(:, 2), myScore(:, 3), 'color', colormap(i, :));
    pause(0.1);
end


end
%
xlabel('PC1');  ylabel('PC2');  zlabel('PC3');
% axis([-2 8 -1 9 -4 4]);
view([121 20]);
hold off


% save figure
savefig(gcf, fullfile('results', 'figures', 'PCA_dynamics_FigureS3D.fig'));
