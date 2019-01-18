addpath(genpath(pwd));

% set the info of data
% global path
path = fullfile('./', 'data', './');  
neuron  = {'Or42a'};
title_str =  'Or42a-10-73pentanol';

% set the parameters about the filters
lfilter = 937; %30s, at the time resolution of 32ms.
lunstable = 1900; %depend on the output data.
rend = 6560; %end of the sequence tploto be analyzed. 
time_point = 0.03203; %time resolution of the recording, unit: s.

% load the stimulus info, including 'ndye' and 'image_times_dye'
load(fullfile('./', 'data', 'Calibration_1023mseq_v1.mat'));

%% Figure 5A, plot one output seq & the input seq, calculat the filter

fname = fullfile(path, 'Or42a_3pentanol_mseq_1.mat'); 
load(fname);
output = normalized_signal{1};
output_n = output/max(output); %normalize the data
fh1 = figure(1);
set(fh1, 'Position', [75, 690, 1200, 300], 'color', 'white');
set(gca,'FontName','Arial'); set(gca,'FontSize',12);
plot(image_times, output_n, 'k'); xlabel('Time(s)');  hold on;
plot(image_times_dye, ndye*0.1 -0.1, 'r'); axis tight; legend('Response', 'Stimulus');  hold off; 
title([title_str, '- Input Output']);


%% Figure S5A part1, Plot the fluorescent traces side by side, check robustness, average all the data
disp('Part 1, subplot all the curves.');

count =0; 

listing = dir('./data/Or42a*.mat');
dffPool = [];

for ii=1:length(listing)
    filename = listing(ii).name;
    fname = [path filename]; 
    load(fname, 'normalized_signal', 'image_times');
    for i =1:length(normalized_signal)
        count = count +1;
        output = normalized_signal{i}; %pick the output signal
        output_s = output(lunstable: rend); %cut off the signal
         
        dffPool(count, :) = output_s'; %pick the output signal

    end
end

fh2 = figure(2);
set(fh2, 'Position', [100, 360, 600, 600], 'color', 'white');
set(gca,'FontName','Arial'); set(gca,'FontSize',12); hold on;

dffRange(:, 1) = min(dffPool, [], 2);
dffRange(:, 2) = max(dffPool, [], 2);

dffStep = dffRange(:, 2) - dffRange(:, 1);
dffGrid = cumsum(dffStep(end: -1: 1));
dffGrid = [0; dffGrid];

for i = 1:7
    dffCurrent = dffPool(i, :)+dffGrid(8-i);    
    plot(image_times_dye(lunstable: rend), dffCurrent, 'k'); hold on;
end
axis tight
xlabel('Time(s)'); 


%% Figure 5B, calculate the linear filters and Figure S5A part 2, filter for each animal
filter = [];
clear filterCut

ndye_n = (ndye - min(ndye))/(max(ndye)-min(ndye));
ndye_s_n = ndye_n(lunstable : rend);

filterWin = 11;

ftimeStep = 1 : 500+125+1;
ftimeStep = ftimeStep - 126 ; 
ftimeStep = ftimeStep + (filterWin-1)/2;  %account for the delay from the filter
ftime = time_point*ftimeStep;

[~,idxL] = min(abs(ftime-(-3)));
[~,idxR] = min(abs(ftime-15));

ftime = ftime(idxL:idxR);

for i = 1 : 7
    filter(:,i) = getfilter2(ndye_s_n, dffPool(i, :)', 125, 500);
    filter_s(:,i) = smooth(filter(:,i), filterWin);

    filterCut(:,i) = filter_s(idxL:idxR, i); 
end

fRange(1, :) = min(filterCut, [], 1);
fRange(2, :) = max(filterCut, [], 1);

fStep = fRange(2, :) - fRange(1, :);
fGrid = cumsum(fStep(end: -1: 1));
fGrid = [0  fGrid];

fh3 = figure(3);
set(fh3, 'Position', [100, 360, 600, 600], 'color', 'white');
set(gca,'FontName','Arial'); set(gca,'FontSize',12); hold on;

for i = 1:7
    fCurrent = filterCut(:, i) + fGrid(8-i);    
    plot(ftime, fCurrent, 'k');  
end
axis tight 
ax = gca;
axis([-3 15 ax.YLim]);
xlabel('Time(s)'); xticks(-3:3:15); 

% normalize all the fitlers and plot on top of each other
filterNorm = filterCut./repmat(max(filterCut,[], 1), [size(filterCut,1) 1]);
figure; hold on;

for i = 2:size(filterNorm, 2)
    plot(ftime, filterNorm(:, i), 'Color', [0.5 0.5 0.5]);
end

plot(ftime, filterNorm(:, 1), 'k');

plot([0 0], [-0.25 1], '-');
axis([-3 15 -0.25 1]); xticks(-3:3:15); yticks(-0.25:0.25:1);

%% Part 5, check the linearity

datalen = rend - lunstable + 1;
pred_data = zeros(datalen, 1);


f_r = getfilter2(ndye_s_n, dffPool(1, :)', 0, lfilter);
f_r = f_r(end : -1 : 1, 1);


for i = 1 : datalen
    ii = ndye_n(i+lunstable-lfilter : i+lunstable);
    pred_data(i) = f_r'*ii;
end

%% Figure S5B, fit the scatter plot using a non-linear function
fh7 = figure(7);
set(fh7, 'Position', [225, 430, 420, 420], 'color', 'white'); 
set(gca,'FontName','Arial'); set(gca,'FontSize',12);

scatter(pred_data, dffPool(1, :), '.k');

axis([-0.2 2.5 -0.2 2.5]);
xlabel('Prediction'); ylabel('Measured');
hold off

x = pred_data;
y = dffPool(1, :)';

ft = fittype('objfunFigS5B(x,a,b,c,d)');
fit_result = fit( x, y, ft);

% When the 'fit' function out of work. Use 'cftool' would be better.
% For Or42a, set the parameters here (from cftool)
% AS AN EXAMPLE, CHANGE THE FOLLOWING PARAMETERS FOR FUTURE RUN
% fit_result.a = 1.274;
% fit_result.b = 3.83;
% fit_result.c = 0.3948;
% fit_result.d = -0.2422;

%draw the fit curve on the scatter plot
xx = min([x;y]) : 0.01 : max([x;y]); 
yy = objfunFigS5B(xx, fit_result.a, fit_result.b, fit_result.c, fit_result.d);
hold on; plot(xx, yy,'r'); hold off;

%% Figure S5A, part 3, compare LN's prediction with actual data

dffPred = zeros(7, datalen);
dffPred(1, :) = objfunFigS5B(x , fit_result.a, fit_result.b, fit_result.c, fit_result.d);

for i = 2: 7
    dffPred(i,:) = (dffPred(1, :)-min(dffPred(1, :))) * dffStep(i)/dffStep(1);
end

figure(fh2);
for i = 1:7
    dffCurrent = dffPred(i, :)+dffGrid(8-i);    
    plot(image_times_dye(lunstable: rend), dffCurrent, 'r'); hold on;
end
axis tight
xlabel('Time(s)'); 