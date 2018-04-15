function C_hat = compute_clean_C(Bcap,alpha,lambda)

% Code to compute cleaned correlation matrix
% 
% Last updated: 2017-03-25
% 
% INPUT
%     Bcap                      Input MSA with phylogenetic effects removed
%     alpha                     Number of significant eigenvectors
%     lambda                    Sample eigenvalues of correlation matrix
%                               formed using Bcap
%
% OUTPUT
%     C_hat                     Cleaned and standardized correlation matrix

%%

[~,M] = size(Bcap);

C = corrcoef(Bcap);
Q = eig_sort(C);

zeta = (M-sum(lambda(1:alpha)))/(M-alpha);
lambda_k_star = [lambda(1:alpha) ; ones(M-alpha,1)*zeta];
C_hat_star = Q*diag(lambda_k_star)*(Q).';
std_C = sqrt(diag(C_hat_star));
C_hat = C_hat_star./(std_C*std_C.');