function [PC, k] = CorrITSPCA(X, gamma_k, num_eig, tol)

% CorrITSPCA --- Correlation based Iterative Thresholding Sparse PCA 
% Code to generate sparse principal components of the given data matrix
% 
% (adapted from Zongming Ma's [1] code
% available at 
% http://www-stat.wharton.upenn.edu/~zongming/software/SPCALab/SPCALab.zip)
% [1] Ma, Z. (2013) Sparse principal component analysis and iterative 
% thresholding. Ann Statist, Vol.41, pp.772-801. 
%
% Last updated: 2017-03-25
% 
% INPUT
%     X           N-by-M data matrix
%     gamma_k     Sparsity parameter
%     num_eig     Number of eigenvectors to be computed
%     tol         Tolerance level in determining convergence
%     
% OUTPUT
%     PC          M-by-(num_eig) matrix, columns corresponding to sparse 
%                 estimate of leading eigenvectors
%     k           1-by-(num_eig) vector recording support sizes of PCs

%%

[n,p] = size(X);

if nargin < 4
    tol = 1/(n^2);
end

%% Standardizing the MSA
% After this step, covariance and correlation of MSA are equal

for kk = 1:p
    X(:,kk) = (X(:,kk)-mean(X(:,kk)))/sqrt(var(X(:,kk)));
end

%% Obtaining the principal components and the corresponding eigenvalues

pn = max(p, n);

[wcPC,D] = eig(cov(X));
d = diag(real(D));
d = flipud(d);
d = d(1:num_eig);
% d = max(d, ones(size(d)));


wcPC = fliplr(wcPC);
wcPC = wcPC(:,1:num_eig);

%% Iterative thresholding algorithm

Q_old = wcPC;
Q_cur = Q_old;
nstep = 0;
S = cov(X);
th = gamma_k;
thresh = @hard_thresholding;

dist = 1;

while (nstep == 0 || dist^2 > tol )
    Q_old = Q_cur;
    Q_cur = S * Q_old;
    for j = 1:num_eig
        Q_cur(:,j) = thresh(Q_cur(:,j), th .* d(j));
    end
    % in case we obtain zero-vector!
    if any(find(Q_cur))==0
        Q_cur = Q_old;
        break
    end
    [Q_cur, ~] = qr(Q_cur, 0);
    nstep = nstep + 1;
    dist = subsp_dist(Q_old, Q_cur);
    if dist == 1
        warning('ITSPCA:OverThresh',...
            'Singularity caused by thresholding!');
    end
end

%% Computing the support of each significant principal component

k = zeros(1, num_eig);
for j = 1:num_eig
    k(j) = length(find(Q_cur(:,j)));
end

%% Sparse principal components

PC = Q_cur;
