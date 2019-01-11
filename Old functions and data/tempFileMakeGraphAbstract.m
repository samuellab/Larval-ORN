% run 'plotHeatMap.m', stop it at line 115, wait until i = 5.

data2DTemp= zeros(17, 11); 
index_row = zeros(17, 11);  index_col = zeros(17, 11);

for ii = 1:17
    for jj = 1:11
        rows = (ii-1)*2 + 1 : ii*2 ;    cols = (jj-1)*2 + 1 : min(jj*2, 21);
        A = data2D(rows, cols);
        [~, iii] = max(A(:)); A(iii) = -1;
        [data2D_temp(ii, jj), indexLocal]= max(A(:));
        [I, J] = ind2sub(size(A), indexLocal);
        index_row(ii, jj) = rows(I);    index_col(ii, jj) = cols(J);
    end
end

figure; imagesc(data2D_temp);
set(gcf, 'Position', [2000 10 560 900]); 
colormap(jet); caxis([0 5.5]);

%% use the row and col index to get elements in other concentration levels

% run 'plotHeatMap.m', stop it at line 115, wait until i = 1.
data2DTemp= zeros(17, 11); 

for ii = 1:17
    for jj = 1:11
        data2D_temp(ii, jj) = data2D(index_row(ii, jj), index_col(ii, jj));
    end
end

figure; imagesc(data2D_temp);
set(gcf, 'Position', [2000 10 560 900]); colormap(jet); caxis([0 5.5]);

%% following processes:
% open in illustrator, scale the x-vs-y to make each box is a square. then
% save to tif figure.

%% shear the image
 
cd 'C:\Users\Lab Admin\Desktop'

fileName = {'1.tif'; '2.tif'; '3.tif'};

A = imread(fileName{2});
tform = affine2d([1 -0.5 0; 0 1 0; 0 0 1]);
J = imwarp(A,tform, 'FillValue', 255); figure; imshow(J);

%%  deal with the sensitivity matrix

% cd 'C:\Users\Lab Admin\Documents\GitHub\larval_olfaction' ;

% run 'AnalyzeEC50matrix.m', stop at line 42.
ecmm = zeros(17, 11);

for ii = 1:17
    for jj = 1:11
        ecmm(ii, jj) = ec50Map(index_row(ii, jj), index_col(ii, jj));
    end
end
figure; imagesc(ecmm); 
cmp = colormap(jet); cmp(1,:) = [0 0 0];
colormap(cmp);
set(gcf, 'Position', [2000 10 560 900]);