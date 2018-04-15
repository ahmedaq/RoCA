%%
% *CODE FOR IDENTIFYING CO-EVOLUTIONARY SECTORS IN ANY PROTEIN USING RoCA*

%% Setting up paths (of functions and data files required) and necessary parameters

clear all;close all;clc

% Adding paths of required functions and datafiles

addpath functions
addpath datafiles

% Setting font type and size

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20)


%% Preprocessing the data

%%%%%%%%%%%%Load MSA matrix%%%%%%%%%%%%%%%%%%%%%%
% Load the amino acid MSA matrix named 'msa' here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[B,Bcap,lambda,site_freq_no_mutation,true_indices,freq_bin,prev_aa,N,M,ls] = preprocessing_gag(msa);

%% Computing the number of significant eigenvectors, alpha

N_shuffles = 1e6; %Number of shuffles
thresh = 1; %percentile to define threshold for significant positive and negative correlations
[lambda_max_rnd,pos_thresh,neg_thresh] = ...
    computing_lambda_max_rnd(Bcap,N_shuffles,thresh);
    
alpha = sum(lambda>lambda_max_rnd);

%% Forming sectors using RoCA
    
% Computing the sparsity threshold, gamma_k for Corr-ITSPCA

gamma_k = computing_gamma_k(alpha,lambda,N,M);

% Computing sparse principal components using Corr-ITSPCA

PC_roca = CorrITSPCA(Bcap, gamma_k, alpha);

% Forming sectors

[sec_eig_roca, sec_eig_roca_true, sec_eig_roca_incl_cs, length_sec_roca] = ...
    form_sectors_roca(PC_roca,alpha,site_freq_no_mutation,true_indices,ls);
    
n_secs_roca = length(sec_eig_roca);
    

%% Computing correlation matrices

C = corrcoef(Bcap); %sample correlation matrix
C_hat = compute_clean_C(Bcap,alpha,lambda); %cleaned standardized correlation matrix

%% Heat maps

heatmap_corr_matrix(C,'eastoutside','Sample correlation matrix')

% Heat map of cleaned correlation matrix 
per_sites_in_roca_sectors = heatmap_corr_matrix_cleaned(C_hat,sec_eig_roca,sec_eig_roca_incl_cs,n_secs_roca,ls,'RoCA');

%% Statistics of RoCA sectors

[mc,mean_abs_corr,per_neg_corr,per_pos_corr] = ...
    stats_sectors(C_hat,sec_eig_roca,n_secs_roca,freq_bin,pos_thresh,neg_thresh);

%% Biplots of eigenvectors

markersize = 5; %Size of circles
jitter_roca = 0.005; %Random jitter to avoid super-imposition of data points

generate_biplots(PC_roca,sec_eig_roca,alpha,markersize,jitter_roca,'RoCA');

%% NED calculation plot

[NED_abs_random,NED_abs_inter,NED_abs_intra] = calculate_NED(B,PC_roca,alpha);
figure_NED(alpha,NED_abs_random,NED_abs_inter,NED_abs_intra)


