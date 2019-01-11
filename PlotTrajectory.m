% folderList = {'\\labnas100\guangwei\Data\20180426\', ...
%     '\\labnas100\guangwei\Data\20180405\'};
% 
% subFolder = {'Orco-4p2b_run101.h5_left_aligned', ...
%     'Orco-4p2b_run202.h5_right_aligned', ...
%     'Orco-4p2b_run301.h5_left_aligned',...
%     'Orco-4p2b_run401.h5_right_aligned', ...
%     'Orco-4p2b_run501.h5_left_aligned', ...
%     'Orco-4p2b_run702.h5_left_aligned', ...
%     'Orco_ia_run101.h5_left_aligned', ...
%     'Orco_ia_run201.h5_left_aligned', ...
%     'Orco_ia_run301.h5_left_aligned', ...
%     'Orco_ia_run401.h5_left_aligned', ...
%     'Orco_ia_run501.h5_left_aligned', ...
%     'Orco_ia_run701.h5_left_aligned'};

folderList = {'\\labnas100\guangwei\Data\20180322\', ...
    '\\labnas100\guangwei\Data\20180330\'};

% subFolder = {'Orco_tdt_ben_run1001.h5_right_aligned', ...
%     'Orco_tdt_ben_run1101.h5_left_aligned', ...
%     'Orco_tdt_ben_run1201.h5_right_aligned',...
%     'Orco_tdt_ben_run601.h5_right_aligned', ...
%     'Orco_tdt_ben_run701.h5_right_aligned', ...
%     'Orco_tdt_ben_run901.h5_left_aligned', ...
%     'Orco_eb_run101.h5_right_aligned', ...
%     'Orco_eb_run401.h5_right_aligned', ...
%     'Orco_eb_run501.h5_right_aligned', ...
%     'Orco_eb_run601.h5_left_aligned', ...
%     'Orco_eb_run701.h5_right_aligned', ...
%     'Orco_eb_run802.h5_left_aligned'};

subFolder = {'Orco_tdt_ben_run601.h5_right_aligned', ...
    'Orco_eb_run101.h5_right_aligned'};


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
for ff = 1:length(subFolder)
fileName = [subFolder{ff}, '_CNMF.mat'];

if ff<=1
    folder = folderList{1};
    colormap = colormap1;
else
    folder = folderList{2};
    colormap = colormap2;
end
fullFileName = fullfile(folder, subFolder{ff}, fileName);
load(fullFileName);

dff = f.CNMF.C_df;
t = f.t;
id = f.neuronID;
oTInfo = f.odor_seq.tcum;

ps = 2:2:10;
ptStart = oTInfo(ps);
ptEnd = oTInfo(ps+1) + 15;

%% load ORN list
load('.\data\doseResponseData.mat', 'ORNList');
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
    a = a - offsetVec;
    
    myScore = a/(coeff');
    
    scatter3(myScore(:, 1), myScore(:, 2), myScore(:, 3), ...
        'MarkerFaceColor', colormap(i, :), ...
        'MarkerEdgeColor', 'w'); 
    plot3(myScore(:, 1), myScore(:, 2), myScore(:, 3), 'color', colormap(i, :));
    pause(0.1);
end


end
%%
xlabel('PC1');  ylabel('PC2');  zlabel('PC3');
% axis([-2 8 -1 9 -4 4]);
view([1,2,1]);
hold off

%%
subFolder = { ...
%     'Orco_tdt_ben_run1001.h5_right_aligned', ...
%     'Orco_tdt_ben_run1101.h5_left_aligned', ...
%     'Orco_tdt_ben_run1201.h5_right_aligned',...
%     'Orco_tdt_ben_run601.h5_right_aligned', ...
%     'Orco_tdt_ben_run701.h5_right_aligned', ...
%     'Orco_tdt_ben_run901.h5_left_aligned', ...
%     'Orco_eb_run701.h5_right_aligned', ...
    'Orco_eb_run101.h5_right_aligned', ...
    'Orco_eb_run401.h5_right_aligned', ...
    'Orco_eb_run501.h5_right_aligned', ...
    'Orco_eb_run601.h5_left_aligned', ...
    'Orco_eb_run802.h5_left_aligned'};



figure; hold on;

for ff = 1:length(subFolder)
fileName = [subFolder{ff}, '_CNMF.mat'];

if ff<=0
    folder = folderList{1};
    colormap = colormap1;
else
    folder = folderList{2};
    colormap = colormap2;
end
fullFileName = fullfile(folder, subFolder{ff}, fileName);
load(fullFileName);

dff = f.CNMF.C_df;
t = f.t;
id = f.neuronID;
oTInfo = f.odor_seq.tcum;

ps = 2:2:10;
ptStart = oTInfo(ps);
ptEnd = oTInfo(ps+1) + 15;

% load ORN list
load('.\data\doseResponseData.mat', 'ORNList');
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

%
i = 5;
tStart = ptStart(i); tEnd = ptEnd(i);
[~, iS]= min(abs(t-tStart)); [~, iE]= min(abs(t-tEnd));
temp= dff(:, iS-1 : iE);

box = zeros(21, size(temp, 2));
for j = 1:length(lut)
    box(lut(j), :) = temp(j, :);
end

dffSeq{i} = box;

a = dffSeq{i};
a = a';
a = a - offsetVec;

myScore = a/(coeff');

%     scatter3(myScore(:, 1), myScore(:, 2), myScore(:, 3), ...
%         'MarkerFaceColor', colormap(i, :), ...
%         'MarkerEdgeColor', 'w'); 
%     plot3(myScore(:, 1), myScore(:, 2), myScore(:, 3), 'color', colormap(i, :));
    scatter3(myScore(:, 1), myScore(:, 2), myScore(:, 3), ...
    'MarkerFaceColor', colormap(6-ff, :), ...
    'MarkerEdgeColor', 'w'); 
plot3(myScore(:, 1), myScore(:, 2), myScore(:, 3), 'color', colormap(6-ff, :));
pause(0.1);


end
%
xlabel('PC1');  ylabel('PC2');  zlabel('PC3');
% axis([-2 8 -1 9 -4 4]);
view([1,2,1]);
hold off