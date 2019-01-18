function f  = getfilter2(ii, oo, l_minus, l_plus)
%% Calcult the filter according to input, output and filter length
% deconvolution process
% the convolution is defined as: o = i * f;
% for discrete system, can be rewritten as o = I * f.
% o: a vector of output, subscript is from 't_plus' to the end
% I: a matrix of input, size is: 'end - t_plus' by 't_minus + t_plus'
% f: a vector of filter, length is 't_minus + t_plus'. 
% Input:
% i: input stimuli
% o: recorded output
% t_minus: length of the minus part of the filter
% t_plus:  length of the plus  part of the filter
% Output:
% f:  raw filter.
% nf: normalized filter

%% 
% define the vector of output, o
o = oo(l_plus + 1: end-l_minus);
 
% define the matrix of I
lfilter = l_minus + l_plus + 1;
I = zeros(length(ii) - lfilter + 1, lfilter);
for j = 1 : length(ii) - lfilter+1
    I(j,:) = ii(lfilter + j - 1: -1 : j)';
end

% calculate the filter
f = I\o;
end