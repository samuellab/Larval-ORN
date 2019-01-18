% set the info of data
% global path
path = fullfile('./', 'data', './');
% n_files = 6;       
neuron  = {'Or42a'};
title_str =  'Or42a-10-73pentanol';

% set the parameters about the filters
lfilter = 937; %30s, at the time resolution of 32ms.
% lfilter = 468; %15s, in case there is a drop on the filter
lunstable = 1900; %depend on the output data.
rend = 6560; %end of the sequence tploto be analyzed. 
% len = 1861;  %half of the seq, used for comparing.
time_point = 0.03203; %time resolution of the recording, unit: s.

% load the stimulus info, including 'ndye' and 'image_times_dye'
% load('Calibration_4C_1203mseq_v1_mid.mat');
load(fullfile('./', 'data', 'Calibration_1023mseq_v1.mat'));

%% Part 0, plot one output seq & the input seq, calculat the filter
disp('Part 0, plot one sequence with stimulus.');
disp('Select one *.mat data to plot.');
% [filename,pathname]  = uigetfile({[path, '*.mat']});
% fname = [pathname filename]; 
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


%% Part 1, Plot the fluorescent traces side by side, check robustness, average all the data
disp('Part 1, subplot all the curves.');

count =0; 
% sum_data = zeros(rend-lunstable+1, 1);

listing = dir('./mseqData/Or42a*.mat');
dffPool = [];

for ii=1:length(listing)
%     [filename,pathname]  = uigetfile({[path, '*.mat']}); 
    filename = listing(ii).name;
    fname = [path filename]; 
    load(fname, 'normalized_signal', 'image_times');
    for i =1:length(normalized_signal)
        count = count +1;
        output = normalized_signal{i}; %pick the output signal
        output_s = output(lunstable: rend); %cut off the signal
         
        dffPool(count, :) = output_s'; %pick the output signal
        
%         output = normalized_signal{i}; %pick the output signal
%         output_re = myRescale(image_times_dye, output, image_times); %resacle the output singal
%         output_re_s = output_re(lunstable: rend); %cut off the signal
%         output_re_s_n = output_re_s/max(output_re_s); %normalize signal
%         
%         sum_data = sum_data + output_re_s_n; %add new data to the sum


%         if count > n_files
%             temp = count;
%         else 
%             temp = n_files;
%         end
%         subplot(7,1,count);
%         plot(image_times_dye(lunstable: rend), output_re_s_n, 'k'); axis tight;
%         if count ==1
%             title([title_str, '- Curve list']);
%         end
    end
end
% set(fh2, 'name', [title_str, '- Curve list'], 'NumberTitle','off');
% output_ave = sum_data/count;
% output_ave_n = output_ave/max(output_ave);

% fh3 = figure(3); %show the averaged normalized output seq.
% set(fh3, 'Position', [125, 640, 1200, 300], 'color', 'white'); 
% set(gca,'FontName','Arial'); set(gca,'FontSize',12);
% plot(image_times_dye(lunstable: rend), output_ave_n, 'k'); 
% axis tight; xlabel('Time(s)'); title([title_str, '- Averaged output curves']);

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


% %% Part 2, Compare filters from different parts
% f  = zeros(lfilter, 1); %filter
% nf = zeros(lfilter, 1); %normalized filter
% 
% % Case 1, whole (from the averaged output seq, the same for the followings)
% ndye_s = ndye(lunstable : rend);
% ndye_s_n = (ndye_s - min(ndye_s))/(max(ndye_s)-min(ndye_s));
% [f, nf] = getfilter(ndye_s_n, output_ave_n, lfilter, 0);
% nf_sm = smooth(nf); nf_sm_n = nf_sm/max(nf_sm);
% 
% fh4 = figure(4);     
% set(fh4, 'Position', [150, 495, 560, 420], 'color', 'white'); 
% set(gca,'FontName','Arial'); set(gca,'FontSize',12);
% plot((0:lfilter-1) * time_point, nf_sm_n, 'r'); hold on; xlabel('lag(s)');
% 
% % Case 2, 1st half
% [f, nf] = getfilter(ndye_s_n(1:lfilter+len), output_ave_n(1:lfilter+len), lfilter, 0);
% nf_sm = smooth(nf); nf_sm_n = nf_sm/max(nf_sm);
% plot((0:lfilter-1) * time_point, nf_sm_n, 'g');
% 
% % Case 3, 2nd half
% [f, nf] = getfilter(ndye_s_n(len+1:end), output_ave_n( len+1:end), lfilter, 0);
% nf_sm = smooth(nf); nf_sm_n = nf_sm/max(nf_sm);
% plot((0:lfilter-1) * time_point, nf_sm_n, 'b'); 
% 
% title([title_str, '- Compare the filters from didifferent parts']);
% legend('Whole', '1st half', '2nd half');
% hold off;

% %% Part 3: Apply the filter from the 1st half of the data to the 2nd half of the data
% [f, nf] = getfilter(ndye_s_n(1:lfilter+len), output_ave_n(1:lfilter+len), lfilter, 0);
% % figure; plot((0:lfilter-1) * time_point, f, 'k'); %show the filter
% 
% new_data = zeros(len,1);
% for i =1:len
%     filter_r = f(end:-1:1);
%     ii = ndye_s_n(len+i+1 : len+i+lfilter);
%     new_data(i) = filter_r'*ii;
% end
% fh5 = figure(5);
% set(fh5, 'Position', [175, 595, 1200, 300], 'color', 'white'); 
% set(gca,'FontName','Arial'); set(gca,'FontSize',12);
% plot(image_times_dye(len+lfilter+1 : len+lfilter+len),...
%     output_ave_n( len+lfilter+1 : len+lfilter+len),'k');
% hold on;
% plot(image_times_dye(len+lfilter+1 : len+lfilter+len),...
%     new_data, 'b'); 
% title([title_str, '- Apply kernel from 1st part to 2nd part']);
% legend('Real', 'Predicted', 'Location', 'NorthEastOutside'); 
% xlabel('Time(s)'); axis tight; hold off;

% %% Part 4: Apply the filter from the 2st half of the data to the 1st half of the data
% [f, nf] = getfilter(ndye_s_n(len+1:end), output_ave_n( len+1:end), lfilter, 0);
% 
% new_data_v2 = zeros(len,1);
% for i =1:len
%     filter_r = f(end:-1:1);
%     ii = ndye_s_n(i+1 : i+lfilter);
%     new_data_v2(i) = filter_r'*ii;
% end
% fh6 = figure(6);
% set(fh6, 'Position', [200, 570, 1200, 300], 'color', 'white'); 
% set(gca,'FontName','Arial'); set(gca,'FontSize',12);
% plot(image_times_dye(lfilter+1 : lfilter+len),...
%     output_ave_n( lfilter+1 : lfilter+len),'k');
% hold on;
% plot(image_times_dye(lfilter+1 : lfilter+len),...
%     new_data_v2, 'b'); 
% legend('Real', 'Predicted', 'Location', 'NorthEastOutside');
% title([title_str, '- Apply kernel from 2nd part to 1st part']);
% xlabel('Time(s)'); axis tight; hold off;


%%
% filter  = zeros(lfilter, 7); %filter
% nf = zeros(lfilter, 7); %normalized filter
filter = [];
clear filterCut
% ndye_s = ndye(lunstable : rend);
% ndye_s_n = (ndye_s - min(ndye_s))/(max(ndye_s)-min(ndye_s));

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
%     nf_sm_n = nf_sm/max(nf_sm);
    
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
%     plot(time_point*(1:size(filterCut, 1)), fCurrent, 'k');  
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

% filterAve = mean(filterNorm, 2);
% plot(ftime, filterNorm(:, i), 'k');
plot([0 0], [-0.25 1], '-');
axis([-3 15 -0.25 1]); xticks(-3:3:15); yticks(-0.25:0.25:1);

%% Part 5, check the linearity
% [f, nf] = getfilter(ndye_s_n, output_ave_n, lfilter, 0);

% new_data_v3 = zeros(2*len,1);
% for i =1:2*len
%     filter_r = f(end:-1:1);
%     ii = ndye_s_n(i+1 : i+lfilter);
%     new_data_v3(i) = filter_r'*ii;
% end

datalen = rend - lunstable + 1;
pred_data = zeros(datalen, 1);


f_r = getfilter2(ndye_s_n, dffPool(1, :)', 0, lfilter);
f_r = f_r(end : -1 : 1, 1);

% idx0 = find(ftime==0);
% countL = idx0 -1;
% countR = length(filter_r) - idx0;

% lunstable = 1900; %depend on the output data.
% rend = 6560; %end of the sequence tploto be analyzed. 

for i = 1 : datalen
    ii = ndye_n(i+lunstable-lfilter : i+lunstable);
    pred_data(i) = f_r'*ii;
end

%% plot all the output data
fh7 = figure(7);
set(fh7, 'Position', [225, 430, 420, 420], 'color', 'white'); 
set(gca,'FontName','Arial'); set(gca,'FontSize',12);

% for i = 1 : datalen
%     plot(pred_data(i), dffPool(1, :), '.'); hold on
% end

scatter(pred_data, dffPool(1, :), '.k');


% axis([-0.1 1.2 -0.1 1.2]);
% plot(-0.1:0.1:1.2, -0.1:0.1:1.2, 'k');
axis([-0.2 2.5 -0.2 2.5]);
xlabel('Prediction'); ylabel('Measured');
% title([title_str, '- Real vs Predicted']); 
hold off

% figure; 
% plot(1:length(pred_data), pred_data, 'k'); hold on;
% plot(1:length(pred_data), dffPool(1, :), 'r'); hold off


%% Part 6, fit the scatter plot using a non-linear function
% x = new_data_v3(1:2*len);
% y = output_ave_n(lfilter+1:lfilter+2*len);
x = pred_data;
y = dffPool(1, :)';

ft = fittype('objfun(x,a,b,c,d)');
fit_result = fit( x, y, ft);

% When the 'fit' function out of work. Use 'cftool' would be better.
% For Or42a, set the parameters here (from cftool)
% JUST AN EXAMPLE, CHANGE THE FOLLOWING PARAMETERS FOR FUTURE RUN
% fit_result.a = 1.274;
% fit_result.b = 3.83;
% fit_result.c = 0.3948;
% fit_result.d = -0.2422;

%draw the fit curve on the scatter plot
xx = min([x;y]) : 0.01 : max([x;y]); 
yy = objfun(xx, fit_result.a, fit_result.b, fit_result.c, fit_result.d);
hold on; plot(xx, yy,'r'); hold off;

%% compare LN's prediction with actual data
% amp_f = max(filterCut, [], 1);

dffPred = zeros(7, datalen);
dffPred(1, :) = objfun(x , fit_result.a, fit_result.b, fit_result.c, fit_result.d);

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


% %% Part 7, apply LN-NLN, compare with real data.
% [f, nf] = getfilter(ndye_s_n(1:lfilter+len), output_ave_n(1:lfilter+len), lfilter, 0);
% 
% new_data = zeros(len,1);
% for i =1:len
%     filter_r = f(end:-1:1);
%     ii = ndye_s_n(len+i+1 : len+i+lfilter);
%     temp = filter_r'*ii;
%     new_data(i) = objfun(temp, fit_result.a, fit_result.b, fit_result.c, fit_result.d);
%     
% end
% figure(fh5); hold on;
% plot(image_times_dye(len+lfilter+1 : len+lfilter+len),...
%     new_data, 'r'); 
% legend('Real', 'Predicted by LN', 'Predicted by LN-NL', 'Location', 'NorthEastOutside');  hold off;

% %% Part 8,  compare fitlers from different trails
% count = 0;
% dataPool = zeros(rend-lunstable+1, 7);
% 
% for ii=1:n_files
%     [filename,pathname]  = uigetfile({[path, '*.mat']}); fname = [pathname filename]; load(fname);
%     for i =1:length(normalized_signal)
%         output = normalized_signal{i}; %pick the output signal
%         output_re = rescale(image_times_dye, output, image_times); %resacle the output singal
%         output_re_s = output_re(lunstable: rend); %cut off the signal
% %         output_re_s_n = output_re_s/max(output_re_s); %normalize signal
% 
%         count = count +1;
%         dataPool(:, count) =  output_re_s; %add new data to the sum
%     end
% end

% %%
% dffRange(1, :) = min(dataPool, [], 1);
% dffRange(2, :) = max(dataPool, [], 1);
% 
% dffStep = dffRange(2, :) - dffRange(1, :);
% dffGrid = cumsum(dffStep(end: -1: 1));
% dffGrid = [0  dffGrid];
% 
% figure; 
% for i = 1:7
%     dffCurrent = dataPool(:, i)+dffGrid(8-i);    
%     plot(image_times_dye(lunstable: rend), dffCurrent, 'k'); hold on;
% end


% %%
% % filter  = zeros(lfilter, 7); %filter
% % nf = zeros(lfilter, 7); %normalized filter
% 
% clear filterCut
% ndye_s = ndye(lunstable : rend);
% ndye_s_n = (ndye_s - min(ndye_s))/(max(ndye_s)-min(ndye_s));
% 
% for i = 1 : 7
%     filter(:,i) = getfilter2(ndye_s_n, dataPool(:, i), 125, 500);
%     filter(:,i) = smooth(filter(:,i));
%     filterCut(:,i) = filter(10:610,i);
% %     nf_sm_n = nf_sm/max(nf_sm);
%     
% end
% 
% %%
% fRange(1, :) = min(filterCut, [], 1);
% fRange(2, :) = max(filterCut, [], 1);
% 
% fStep = fRange(2, :) - fRange(1, :);
% fGrid = cumsum(fStep(end: -1: 1));
% fGrid = [0  fGrid];
% 
% figure; 
% for i = 1:7
%     fCurrent = filterCut(:, i) + fGrid(8-i);    
%     plot(time_point*(1:size(filterCut, 1)), fCurrent, 'k'); hold on;
% end
% 
% 


%% apply the filter from the 1st trail to other trials, with considering the amplitude differneces

