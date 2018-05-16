load('signal_list.mat');

% find unique odors
[C,ia,ic] = unique(T_rows(:, 1)); 
for i = 1:length(ia)
    row_per_odor = find(ic == i);
    figure;
    imagesc(sigMat(row_per_odor, :)); title(string(C.Odor(i)));
end