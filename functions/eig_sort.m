function [eigvect,lambda] = eig_sort(C)

% Code to generate eigenvalues and eigenvectors arranged in descending 
% order
%
% Last updated: 2017-03-25
% 
% INPUT
%     C                         M-by-M correlation matrix
%
% OUTPUT
%     eigvect                   M-by-M matrix of eigenvectors arranged in
%                               descending order
%     lambda                    M-by-1 vector of eigenvalues arranged in
%                               descending order

%%

[eigvect_unsorted,lambda_unsorted]=eig(C);
[lambda,lambda_order]=sort(diag(lambda_unsorted),'descend');
eigvect=eigvect_unsorted(:,lambda_order);
