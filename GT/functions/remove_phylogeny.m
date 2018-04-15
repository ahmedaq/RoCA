function [msa_without_phylogeny,lambda] = remove_phylogeny(msa_bin,method_remove_phylogeny)

% Code to remove the phylogenetic effects from the input MSA
% 
% Last updated: 2017-03-25
% 
% INPUT
%     msa_bin                   N-by-M binarized MSA
%     method_remove_phylogeny   Method to remove phylogenetic effects
%                               - 0 (linear regression)
%                               - 1 (logistic regression)
%                               - 2 (clipping)
%     
% OUTPUT
%     msa_without_phylogeny     MSA with effect of phylogeny potentially 
%                               removed
%     lambda                    Eigenvalues of the MSA without phylogeny

%%
[N,M] = size(msa_bin);
C = corrcoef(msa_bin);
[Q,lambdaB] = eig_sort(C);
q_1 = Q(:,1);
proj_1 = msa_bin*q_1;
  
if method_remove_phylogeny == 0 
  
    x = zeros(M,2);
    
    parfor k = 1:M
        A = [ones(N,1) proj_1];
        b = msa_bin(:,k);
        x(k,:) = inv(A'*A)*A'*b;        
    end
    
    msa_without_phylogeny = msa_bin - proj_1*x(:,2).';
    C_without_phylogeny = corrcoef(msa_without_phylogeny);

elseif method_remove_phylogeny == 1 
    
    parfor k = 1:M
        b = msa_bin(:,k);
        probs_01 = mnrval(mnrfit(proj_1,b+1),proj_1); 
        %because mnrfit requires greater than 1 in second
        %argument --- two categories 1 and 2
        X(:,k) = probs_01(:,2);
    end
    
    msa_without_phylogeny = msa_bin - X;
    C_without_phylogeny = corrcoef(msa_without_phylogeny);
    
elseif method_remove_phylogeny == 2 
    
    
    lambda_wo_phylogeny = lambdaB;
    lambda_wo_phylogeny(1) = [];
    eigvect_wo_phylogeny = Q;
    eigvect_wo_phylogeny(:,1) = [];
    
    C_without_phylogeny = eigvect_wo_phylogeny*diag(lambda_wo_phylogeny)*eigvect_wo_phylogeny';
    
end

[~,lambda] = eig_sort(C_without_phylogeny);