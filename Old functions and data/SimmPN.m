%% select odor index
totalOdorNum = 18;

% tranMAll = zeros(5, 40, totalOdorNum);
aveFRAll = zeros(5, 40, totalOdorNum);
%% define data file
conFile = fullfile('data', 'ORN-PickyLN-mPN_Left_Merged.xlsx');
matFile = fullfile('data', 'AveRawDataMatrix.mat');
fitFile = fullfile('AnalysisResults', 'FitDoseResponse.mat');

%% load data
addpath(genpath(pwd));
[ connect, synType, inputs, nList ] = GetCon(conFile);
load(matFile, 'dataRawAve', 'infoORNList', 'odorList', 'concList');
load(fitFile, 'aMatrix');

%% Normalize data using fitted amplitude
data = dataRawAve;
% normalize data of each ORN-odor pair using fitted r_max
[r, c] = find(~isnan(aMatrix));
for i = 1:length(r)
    data(r(i), c(i), :) = data(r(i), c(i), :)/aMatrix(r(i), c(i));
end

data(5, :, :) = [];   odorList(5) = []; %remove odor 5
data(data<0) = 0; %set the negative values to 0

%% calculate exp to connnect ORN sequance 
exp2con = GetORNSeq(infoORNList, nList(inputs));

%% define dose-response simulaiton params
pDuration = 0.5;
pNum = length(concList);
numSeg = 2*pNum+1;
params.time = numSeg * pDuration;

%% define LIF model params
params.dt = 1e-4;
params.V0 =-52e-3;
params.spikeMax=20e-3;
params.spikeMin=-72e-3;
params.spikeThr=-45e-3;
params.apDuration=2;
params.bkgFire = 10;
params.maxFire = 200;

%%
draw = 0;

%%
for iii = 1:totalOdorNum
    
    odorIndex = iii;

%% defien AP
apRiseTime=round((params.apDuration*10)/2);
apRise=normpdf(-1:1/apRiseTime:0);
apRise=(apRise-min(apRise))/(max(apRise)-min(apRise));
apRise=apRise*(params.spikeMax-params.spikeThr)+params.spikeThr;
apFallTime=round((params.apDuration*9)/2);
apFall=sin(pi()/2:pi()/apFallTime:3*pi()/2);
apFall=(apFall-min(apFall))/(max(apFall)-min(apFall));
apFall=apFall*(params.spikeMax-params.spikeMin)+params.spikeMin-.0001;
AP=[apRise apFall]';
AP(1)=[];
AP(end+1)=AP(end)+0.0001;
AP(end+1)=AP(end)+0.0001;
apLength=length(AP);

%% generate input seq
pRaw = ones(1, (2*pNum+1) )* pDuration;    %s
N = params.time/params.dt;
VinExp = zeros(N, length(exp2con)) + params.V0;
pRaw = pRaw/params.dt;

inRanges = zeros(length(pRaw), 2);
inRanges(1, 1) = 1; inRanges(end, 2) = N;
for i = 2 :length(pRaw)-1
    inRanges(i, 1) = sum(pRaw(1:i-1))+1;
    inRanges(i, 2) = sum(pRaw(1:i));
end
inRanges(1, 2) = inRanges(2, 1) -1; inRanges(end, 1) = inRanges(end-1, 2)+1;

%
inRaw = data(odorIndex,:,:);
inRaw = reshape(inRaw, [], pNum);
inRawF = inRaw * params.maxFire + params.bkgFire;
inRawFLong = zeros(length(inRawF(:, 1)), numSeg) + params.bkgFire;
inRawFLong(:, 2:2:2*pNum) = inRawF;

for n = 1 : length(exp2con)
    for i = 1 : numSeg
        numSpike = round(pRaw(i) * params.dt * inRawFLong(n, i));
        range = inRanges(i, 2) - inRanges(i, 1)+1;

%         if inRawFLong(n, i) < 1.5 * params.bkgFire
%             temp = round(range/apLength);
%             r = randperm(temp - 1);
%             rP = r(1 : numSpike);
%             for j = 1 : numSpike
%                 startP = rP(j) * apLength + inRanges(i, 1);
%                 VinExp(startP : startP + apLength-1, n) = AP;
%             end
%         else
            interval = round(range/numSpike);
            for j = 1 : numSpike
                r = randperm(interval - apLength);
                startP = interval * (j-1) + r + inRanges(i, 1);
%                 startP = interval * j - apLength + inRanges(i, 1);
                VinExp(startP : startP + apLength-1, n) = AP;
            end
%         end
    end
end

% convert the input matrix to the order in the connect 
inputVs = zeros(N, length(inputs)) + params.V0;
for i = 1:length(exp2con)
    inputVs(:, exp2con(i)) = VinExp(:, i);
end

% find the missed ORNs in the experiment data, add random firing
temp = ismember(1:length(inputs), exp2con);
exIndex = find(temp==0);

VinEx = zeros(N, length(exIndex)) + params.V0;
inRawFLongEx = zeros(length(exIndex), numSeg) + params.bkgFire;

for n = 1 : length(exIndex)
    for i = 1 : numSeg
        numSpike = round(pRaw(i) * params.dt * inRawFLongEx(n, i));
        range = inRanges(i, 2) - inRanges(i, 1)+1;
            interval = round(range/numSpike);
            for j = 1 : numSpike
                r = randperm(interval - apLength);
                startP = interval * (j-1) + r + inRanges(i, 1);
                VinEx(startP : startP + apLength-1, n) = AP;
            end
    end
end

for i = 1:length(exIndex)
    inputVs(:, exIndex(i)) = VinEx(:, i);
end

%%
tic
out = flyLIF(params,connect,synType,inputs,inputVs);
toc

%% plot
if draw == 1
    figure; set(gcf,'units','normalized','position',[ 0.0365 0.0102 1.9589 0.9120]);
    nNum = length(synType);
    colN = 2;
    rowN = nNum/colN;
    for i = 1:length(synType)
        subplot(rowN, colN, i);
        plot((1:length(out.V(1:N,i))).*out.P.dt, out.V(1:N,i)); 
        axis([0 out.P.time -0.08 0.02]); title(nList{i});
    end
%     saveas(gcf, [odorList{odorIndex}, '.fig']);
    close();
end

%% plot the spike cont plot 
tBin = 0.1; %s
sc = v2sc( out.V, out.P.dt, tBin, out.P.spikeMax );

if draw == 1
    figure; set(gcf, 'position', [70   375   1050   600]);
    a = axes;
    imagesc(sc'); 
    set(gca,'XTick', (pDuration/tBin : pDuration/tBin : params.time/tBin) -  pDuration/tBin/2);
    xLStr = cell(1, numSeg);
    xLStr(2:2:(numSeg-1)) = cellstr(num2str(concList));
    set(gca,'XTickLabel',xLStr);
    set(gca,'xaxisLocation','top');
    set(gca,'YTick',1:length(sc(1,:)));
    set(gca,'YTickLabel',nList);
    title(odorList{odorIndex});
    colormap(jet); colorbar;
%     savefig([ odorList{odorIndex}, '_heatMap.fig']);
end
        
%% simplify the data, get transient True/False matrix and ave Firing Rate 
[ ~, aveFR ] = GetFiringRate( sc, tBin, pDuration );
% tranM = tranM(2:2:numSeg-1, :);
aveFR = aveFR(2:2:numSeg-1, :) - aveFR(1:2:numSeg-2, :);

% tranMAll(:, :, iii) = tranM;
aveFRAll(:, :, iii) = aveFR;

 
if draw == 1
%     figure; 
%     set(gcf, 'position', [45   80   350   600]);
% %     a = axes;
%     imagesc(tranM'); 
%     set(gca,'XTick', (1 : 1 : pNum));
%     xLStr= cellstr(num2str(concList));
%     set(gca,'XTickLabel',xLStr);
%     set(gca,'xaxisLocation','top');
%     set(gca,'YTick',1:length(tranM(1,:)));
%     set(gca,'YTickLabel',nList);
%     title(odorList{odorIndex});
%     colormap(jet);
%     savefig([ odorList{odorIndex}, '_TransientCheck.fig']);
    
    figure; 
    set(gcf, 'position', [395   80   400   600]);
    imagesc(aveFR'); 
    set(gca,'XTick', 1 : 1 : pNum);
    xLStr= cellstr(num2str(concList));
    set(gca,'XTickLabel',xLStr);
    set(gca,'xaxisLocation','top');
    set(gca,'YTick',1:length(aveFR(1,:)));
    set(gca,'YTickLabel',nList);
    title(odorList{odorIndex});
    colormap(jet); colorbar;
%     savefig([ odorList{odorIndex}, '_AveFR.fig']);
end
end

%%
ORNAct = aveFRAll(:, 1:21, :);
mPNpiLNAct = aveFRAll(:, 22:end, :);
piLNAct = aveFRAll(:, 22:26, :);
mPNAct = aveFRAll(:, 27:end, :);

mPNAct = mPNAct/max(mPNAct(:));

%% plot mPN activity
odorOrder = [17 12 15 2 10 3 4 16 9 1 18 8 6 5 7 13 14 11]; % the order is consistant to figure 2
mPNList = nList(27:end);

mPN_DataTemp = mPNAct;
[z, neurons, odors] = size(mPN_DataTemp);

mPNnewM = mPN_DataTemp;
for i = 1:odors
    mPNnewM(:,:,i) = mPN_DataTemp(:,:,odorOrder(i));
end

for i = 1 : z
    data2D = mPNnewM(i,:,:);
    data2D = reshape(data2D, [neurons, odors]);
    data2D = data2D';
    
    figure;
    a = axes;
    set(a, 'CLim', [0 1]);
    imagesc(data2D); 
    set(gca,'XTick',1:length(mPNList));
    set(gca,'XTickLabel',mPNList(:));
    set(gca,'xaxisLocation','top');
    set(gca,'YTick',1:length(odorList));
    set(gca,'YTickLabel',odorList(odorOrder));
    ax = gca; ax.XTickLabelRotation = 45;
    title(concList(i));
    colormap(jet);
    caxis([0 1]);
    set(gcf, 'Position', [680   558   510   375]);
%     saveas(gcf, ['C:\Users\Lab Admin\Desktop\mPN', num2str(i), '.svg'] );
end


%% PCA on the mPNs' activity
ORNActTemp = permute(mPNAct,[3 2 1]); 
dataTall = permute(ORNActTemp, [1 3 2]); 
dataTall = reshape(dataTall,[],size(ORNActTemp,2),1);
[coeff,score,latent,tsquared,explained,mu] = pca(dataTall, 'Centered', true);

figure; 
plot(1:length(explained), explained, 'ok');
xlabel('PC'); ylabel('Variance Explained(%)');
axis([0 length(explained) 0 80]);

%bar plot
figure;
bar(1:5, explained(1:5), 0.5, 'k' )
axis([.5 5.5 0 70]);
set(gca,'XTick',1:5);
set(gca,'YTick',0:10:70);
set(gca,'YTickLabel', 0:10:70);
xlabel('Principal Components');
ylabel('% of variance explained');

%% plot the data projection onto the first 3 PCs
%define colormap for odors
odorColorMapRaw = [...
    138 198 64;...	%1-pentanol, 35a
    175 52 147 ;...	%3-pentanol
    0 145 87 ;...	%6-methyl-5-hepten-2-ol, 13a
    78 71 157;...	%3-octanol, 85c
    251 106 74;...	%methyl phenyl sulfide, 24a
    239 59 44; ...	%anisole, 30a
    203 24 29 ; ...	%2-acetylpyridine, 22c
    217,95,14; ...  %2,5-dimethylpyrazine, 33b
    67 182 73; ...  %pentyl acetate, 47a
    118,42,131  ;...%geranyl acetate, 82a
    140,81,10 ; ... %2-methoxyphenyl acetate, 94a
    45 85 166; ...  %trans,trans-2,4-nonadienal, 
    103 0 13;...    %4-methyl-5-vinylthiazole, 45b
    165 15 21; ...  %4,5-dimethylthiazole, 59a
    64,0,75; ...    %4-hexen-3-one, 42a
    30 142 205; ...	%2-nonanone, 45a
    153,112,171; ...%acetal, 42b
    254 192 15; ... %2-phenyl ethanol
    ];
odorColorMap = odorColorMapRaw./255.0;
markerSize = 400;     mksize_temp = [ 0.1 0.2 0.4 0.6 1];

figure;
set(gcf, 'Position', [700 80 1100 850]);
for i =1:length(odorList)
	seqTemp = length(odorList) * [0 1 2 3 4] + i;
	scatter(score(seqTemp, 1), score(seqTemp, 2), ...
         markerSize*mksize_temp, odorColorMap(i,:),'fill'); 
    hold on;
    
	plot(score(seqTemp, 1), score(seqTemp, 2), ...
        'color', odorColorMap(i,:), ...
        'LineWidth', 1.5);
	pause(0.1);
end
axis tight

%add odor name
% dx = 0.1; dy = 0.1;  % displacement so the text does not overlay the data points
% maxInd = (1:length(odorList)) + 4*length(odorList);
% text(score(maxInd, 1) + dx, score(maxInd, 2) +dy, odorList);  
xlabel(['PC', num2str(1),' (', num2str(explained(1), 3), '% of variance)']);
ylabel(['PC', num2str(2),' (', num2str(explained(2), 3), '% of variance)']);