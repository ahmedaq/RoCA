function binaryvect = get_binary_vecs(i)

% Code to generate a matrix of all possible binary vectors given the number
% of bits i
% 
% Last updated: 2017-03-25
% 
% INPUT
%     i                     Number of input bits
%
% OUTPUT
%     binaryvect            Matrix of all possible binary i bit vectors

D = [0:2^i - 1];
binaryvect = double(dec2bin(D))-48;