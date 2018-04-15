function vec_triu = obtain_vec_triu_matrix(C)

% Code to obtain a vector consisting of elements in the upper triangilar 
% part of a matrix
% 
% Last updated: 2017-03-25
% 
% INPUT
%     C                         Input matrix
%
% OUTPUT
%     vec_triu                  vector consisting of elements in the upper
%                               triangilar part of a matrix

%%
vec_triu = [];
for kk = 1:length(C)-1
    for mm = kk+1:length(C)
       vec_triu = [vec_triu C(kk,mm)];
    end
end