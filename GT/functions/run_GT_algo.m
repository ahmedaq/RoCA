function [mean_tpr_spca,mean_fdr_spca,PM_max_spca,...
    mean_tpr_pca,mean_fdr_pca,PM_max_pca] = run_GT_algo(rnm,ITER,msa_bin,r,M,Nbar,supp_actual,verbose)

% Main code running the ground truth simulation on generated binary data
% 
% Last updated: 2017-03-25

%% 

num_eig = r;
msa_orig = msa_bin;
N = M*rnm;

mean_tpr_spca = zeros(ITER,1);
mean_fdr_spca = zeros(ITER,1);
PM_max_spca = zeros(ITER,1);
mean_tpr_pca = zeros(ITER,1);
mean_fdr_pca = zeros(ITER,1);
PM_max_pca = zeros(ITER,1);

for iter = 1:ITER
    
    B = msa_orig(randperm(Nbar,N),:);
    
    freq_single_mutation = mean(B);
    site_freq_no_mutation = find(freq_single_mutation==0);
    
    B(:,site_freq_no_mutation)=[];
    true_indices = setdiff([1:M],site_freq_no_mutation);
       
    %% Removing phylogeny

    [Bcap,lambda] = remove_phylogeny(B,0);
    
    %% Computing the number of significant eigenvectors, alpha
    
    N_shuffles = 20; %Number of shuffles
    thresh = 1; %percentile to define threshold for significant positive and negative correlations
    [lambda_max_rnd,pos_thresh,neg_thresh] = ...
        computing_lambda_max_rnd(Bcap,N_shuffles,thresh);
    
    %Assuming the correct number of eigenvalues is obtained and focus on comparison of PCA and SPCA sectors formed
    alpha = num_eig;
    
    %% Forming sectors using SPCA (Corr-ITSPCA)
    
    % Computing the sparsity threshold, gamma_k
    gamma_k = computing_gamma_k(num_eig,lambda,N,M);
    
    % Computing sparse principal components using SPCA
    PC_spca = CorrITSPCA(Bcap, gamma_k, num_eig);
    
    % Forming sectors using SPCA
    [sec_eig_spca] = ...
        form_sectors_spca(PC_spca,num_eig,site_freq_no_mutation,true_indices,ls);
    
    % Finding the best aligned sectors to units
    [sec_kappa_star_spca] = find_sec_kappa_star(num_eig,sec_eig_spca,supp_actual,M);
    
    % Compute metrics 
    [mean_tpr_spca(iter),mean_fdr_spca(iter),PM_max_spca(iter)] = ...
        compute_metrics(sec_kappa_star_spca,supp_actual,num_eig);
    
    %% Forming sectors using PCA
    
    [~, sec_eig_pca] = ...
        form_sectors_pca(Bcap,num_eig,site_freq_no_mutation,true_indices,M,pos_thresh,neg_thresh);
    
    % Finding the best aligned sectors to units
    sec_eig = cell(1,num_eig);
    for kk = 1:length(sec_eig_pca)
        sec_eig{kk} = sec_eig_pca{kk};
    end
    
    [sec_kappa_star_pca] = find_sec_kappa_star(num_eig,sec_eig,supp_actual,M);
    
    % Compute metrics 
    [mean_tpr_pca(iter),mean_fdr_pca(iter),PM_max_pca(iter)] = ...
        compute_metrics(sec_kappa_star_pca,supp_actual,num_eig);
    
    if verbose == 1
            
        fprintf('iter = %d, TPR_SPCA = %.2f, TPR_PCA = %.2f, FDR_SPCA = %.2f, FDR_PCA = %.2f, PMmax_SPCA = %.2f, PMmax_PCA = %.2f\n',...
            iter,mean_tpr_spca(iter),mean_tpr_pca(iter),mean_fdr_spca(iter),mean_fdr_pca(iter),...
            PM_max_spca(iter),PM_max_pca(iter))
        
    end
    
end