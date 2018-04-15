function [dist] = subsp_dist(A, B)

% Zongming Ma's [1] code available at 
% http://www-stat.wharton.upenn.edu/~zongming/software/SPCALab/SPCALab.zip
% [1] Ma, Z. (2013) Sparse principal component analysis and iterative 
% thresholding. Ann Statist, Vol.41, pp.772-801. 

%%
[A, R] = qr(A, 0);
[B, R] = qr(B, 0);

overlap = svd(A' * B);
overlap = overlap(end);
dist = 1 - overlap^2;
