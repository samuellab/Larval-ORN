function exp2con = GetORNSeq(expStr, conStr)

% simplify the ORN name 
conStrShort = cell(length(conStr), 1);
for i = 1 : length(conStr)
    strTemp = conStr{i};
    indexTemp = strfind (strTemp, ' '); indexTemp = indexTemp(1);
    conStrShort{i} = strTemp(1:indexTemp-1);
end

% find the ORN sequence to make them match each other
exp2con = zeros(length(expStr), 1);
for i = 1: length(exp2con)
    strTemp = expStr{i};
    strTemp = strTemp(3:end);
    if strcmp(strTemp, '33b-47a')
        strTemp = '47a';
    elseif strcmp(strTemp, '94a-94b')
        strTemp = '94a';
    end
    indexC = strfind(conStrShort, strTemp);
    exp2con(i) = find(not(cellfun('isempty', indexC)));
end
end