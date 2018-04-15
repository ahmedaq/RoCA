function [lambda_max_rnd,pos_thresh,neg_thresh] = ...
    computing_lambda_max_rnd(Bcap,N_shuffles,thresh)

%Code to compute the maximum eigenvalue obtained in the null case
% 
% Last updated: 2017-03-25
% 
% INPUT
%     Bcap                  N-by-M MSA 
%     N_shuffles            Number of times to repeat shuffling procedure
%     thresh                percentile to define the significant positive
%                           and negative correlations (e.g., thresh = 1
%                           implies, 1st percetile and 99th percentile will
%                           be used to define the significant negative and
%                           positive correlations, respectively.)
%
% OUTPUT
%     lambda_max_rnd        maximum eigenvalue obtained in the null case
%     pos_thresh            threshold defining significnat positive
%                           correlations
%     neg_thresh            threshold defining significnat negative
%                           correlations

%%

[N_seq,N_pos] = size(Bcap);
lambda_rnd=zeros(N_shuffles,1);
C_random = zeros(N_shuffles*N_pos*(N_pos-1)/2,1);

for s=1:N_shuffles   
    msa_bin_rnd=zeros(N_seq,N_pos);
    for pos=1:N_pos 
        perm_seq=randperm(N_seq);
        msa_bin_rnd(:,pos)=Bcap(perm_seq(:),pos);
    end    
    C_rnd = corrcoef(msa_bin_rnd);
    C_random((s-1)*N_pos*(N_pos-1)/2+1:s*N_pos*(N_pos-1)/2,1) = obtain_vec_triu_matrix(C_rnd);
    [~,lambda_r]=eig_sort(C_rnd);
    lambda_rnd(s) = max(lambda_r);
end

lambda_max_rnd = max(lambda_rnd(:));

C_random = [C_random(:)];
pos_thresh = prctile(C_random,100-thresh);
neg_thresh = prctile(C_random,thresh);