clear all;close all;clc

% Comparing the performance of RoCA and PCA in inferring co-evolutionary
% sectors using binary synthetic data
% 
% Last updated: 2017-03-25
% 
% INPUT
% Inputs can be chaged using the model parameters below.
%
% OUTPUT
% This script produces three boxplot figures that compare the performance
% of RoCA and PCA using the following metrics
%       - Mean true positive rate
%       - Mean false discovery rate
%       - Maximum percentage mismatch

%% Setting paths and font

addpath('functions')
addpath('lib')
addpath('interface')

% Setting font type and size

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20)

%% Model parameters

Nbar = 1e6;                 %number of binary samples to generate
M = 500;                    %number of variables or sites
r = 5;                      %number of spikes
S = (0.14:-0.02:0.06)*M;    %support of eigenvectors
ell_0 = 8;                  %eigenvalue for phylogeny
ell_1 = 6;                  %eigenvalue for the top spiked eigenvector (after phylogeny) 
ell_r = 4;                  %eigenvalue for the smallest spiked eigenvector
% Specifying the 1-point probabilities of variables (that define variance)
max_1pt = 0.2; 
min_1pt = 0.1; 
joint_supp = sum(S); % number of non-zero entries in leading eigenvectors (joint support, characterizing joint sparsity)

fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('Model parameters\n')
fprintf('-----------------------------------------------------------------------------------\n')
fprintf(' Num. samples (Nbar): %d  \n Num. sites (M): %d \n',Nbar,M)
fprintf(' Num. spikes (r): %d  \n Joint support: %d  \n Support sizes (Si): %s \n',r,joint_supp,mat2str(S))

[C,mu,supp_actual] = spiked_corr_model(M,r,S,ell_0,ell_1,ell_r,max_1pt,min_1pt);

%% Generating binary synthetic data using the model parameters

msa_bin = generate_bin_data_EPmethod(Nbar,C,mu);

muHat = mean(msa_bin,1);    % estimate mean
CHat = cov(msa_bin);        % estimate covariance    

fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('Statistics of genrated binary data\n')
fprintf('-----------------------------------------------------------------------------------\n')
fprintf(' Target mean X1:  %.3f     Estimated mean X1:  %.3f\n',mu(1),muHat(1))
fprintf(' Target mean X2:  %.3f     Estimated mean X2:  %.3f\n',mu(2),muHat(2))
fprintf(' Target cov X1X2: %.3f     Estimated cov X1X2: %.3f\n',C(1,2),CHat(1,2))

%% Running SPCA and PCA on generated data and calculating the metrics: Mean TPR, Mean FDR, and PM_max for different N/M

ITER = 50;
verbose = 0; %1: display each iteartion output, 0: suppress display

fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('Running PCA and SPCA on GT data for N/M = 2\n')
fprintf('-----------------------------------------------------------------------------------\n')
[mean_tpr_spca_rnm2,mean_fdr_spca_rnm2,PM_max_spca_rnm2,...
    mean_tpr_pca_rnm2,mean_fdr_pca_rnm2,PM_max_pca_rnm2]...
    = run_GT_algo(2,ITER,msa_bin,r,M,Nbar,supp_actual,verbose);

fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('Running PCA and SPCA on GT data for N/M = 4\n')
fprintf('-----------------------------------------------------------------------------------\n')
[mean_tpr_spca_rnm4,mean_fdr_spca_rnm4,PM_max_spca_rnm4,...
    mean_tpr_pca_rnm4,mean_fdr_pca_rnm4,PM_max_pca_rnm4]...
    = run_GT_algo(4,ITER,msa_bin,r,M,Nbar,supp_actual,verbose);

fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('Running PCA and SPCA on GT data for N/M = 8\n')
fprintf('-----------------------------------------------------------------------------------\n')
[mean_tpr_spca_rnm8,mean_fdr_spca_rnm8,PM_max_spca_rnm8,...
    mean_tpr_pca_rnm8,mean_fdr_pca_rnm8,PM_max_pca_rnm8]...
    = run_GT_algo(8,ITER,msa_bin,r,M,Nbar,supp_actual,verbose);

%% Generating figures

% Boxplor of mean TPR

data_tpr = [mean_tpr_pca_rnm2.' mean_tpr_spca_rnm2.' zeros(1,ITER) ...
    mean_tpr_pca_rnm4.' mean_tpr_spca_rnm4.' zeros(1,ITER) ...
    mean_tpr_pca_rnm8.' mean_tpr_spca_rnm8.'];
ylimdata_tpr = [-.05 1.05];
ytickdata_tpr = 0:.2:1;
ylabeldata_tpr = 'Mean TPR';

fig_boxplot(data_tpr,ITER,ylimdata_tpr,ytickdata_tpr,ylabeldata_tpr)

% Boxplor of mean FDR

data_fdr = [mean_fdr_pca_rnm2.' mean_fdr_spca_rnm2.' zeros(1,ITER) ...
    mean_fdr_pca_rnm4.' mean_fdr_spca_rnm4.' zeros(1,ITER) ...
    mean_fdr_pca_rnm8.' mean_fdr_spca_rnm8.'];
ylimdata_fdr = [-.05 1.05];
ytickdata_fdr = 0:.2:1;
ylabeldata_fdr = 'Mean FDR';

fig_boxplot(data_fdr,ITER,ylimdata_fdr,ytickdata_fdr,ylabeldata_fdr)

% Boxplot of maximum PM 

data_PMmax = [PM_max_pca_rnm2.' PM_max_spca_rnm2.' zeros(1,ITER) ...
    PM_max_pca_rnm4.' PM_max_spca_rnm4.' zeros(1,ITER) ...
    PM_max_pca_rnm8.' PM_max_spca_rnm8.'];
ylimdata_PMmax = [-5 105];
ytickdata_PMmax = 0:20:100;
ylabeldata_PMmax = 'PM_{max}';

fig_boxplot(data_PMmax,ITER,ylimdata_PMmax,ytickdata_PMmax,ylabeldata_PMmax)