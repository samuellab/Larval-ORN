load('signal_list.mat');

%% change the order of the ORNs
target_ORN_list = { '33b/47a'; '45a'; '83a'; '35a'; '42a'; '59a'; '1a';	'45b'; '63a'; ...
    '24a'; '67b'; '85c'; '13a'; '30a'; '82a'; '22c'; '42b';	'33a'; '49a'; '74a'; '94a/94b'};

current_ORN_list = string(ORNs);

mat_reorder = []; %save the data in a new order in this variable
for i = 1:length(target_ORN_list)
    mystr = target_ORN_list{i};
    myIndex = find(strcmp( current_ORN_list , mystr));
    mat_reorder(:, i) = sigMat(:, myIndex);
end

%% define new matrix and table head
mat_new = [];
new_table_var_name = {'Odor', 'Exp_ID', 'Concentration'};
table_new = table([], [], [], 'VariableNames', new_table_var_name);

%% reformat the matrix & table
[C,ia,ic] = unique(T_rows(:, 1));  %find unique odors

for i = 1:length(ia) %go through each odor
    row_per_odor = find(ic == i); % find the rows for each odor
    
%     table_pool = T_rows(row_per_odor, :);
%     mat_pool = mat_reorder(row_per_odor, :);

    % organize the experiment ID
    exp_ID_part1 = T_rows(row_per_odor, 3);
    exp_ID_part2 = T_rows(row_per_odor, 4);
    exp_ID = strcat( string(exp_ID_part1.Date), '_', string(exp_ID_part2.AnimalNumber));
    
    % find unique exp_IDs for each odor
    [ID, ia_ID, ic_ID] = unique(exp_ID);
    
    % cluster the data from the same animal
    for j = 1:length(ia_ID)
        index_exp = find (ic_ID == j); % fidn the index
        
        rows_exp = row_per_odor(index_exp);
 
        % append data to the new matrix
        mat_new = [mat_new; mat_reorder(rows_exp, :)];
        
        % append the table
        table_block = table('Size',[length(rows_exp) 3], 'VariableTypes',{'string','string','double'}, 'VariableNames', new_table_var_name);
        for k = 1:length(rows_exp)
            table_block(k, 1) =  T_rows(rows_exp(k), 1);
            table_block{k, 2} =  ID(j);
            table_block{k, 3} =  str2num(string(T_rows{rows_exp(k), 2}));
        end
        
        table_new = [table_new; table_block];
    end
end

%% visualzie the data
for i = length(ia):-1:1 %go through each odor
     row_per_odor = find(ic == i); % find the rows for each odor
     figure; imagesc(mat_new(row_per_odor, :)); 
     xticks(1:21);
	 xticklabels( target_ORN_list ); xtickangle(-45);
     yticks(1:length(row_per_odor));
     yticklabels( num2str(table_new.Concentration(row_per_odor)) );  
     title(string(C.Odor(i)));
end

%% deal with saturations

