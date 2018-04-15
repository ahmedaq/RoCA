%%
% *CODE FOR IDENTIFYING CO-EVOLUTIONARY SECTORS IN HIV NEF USING SCA*

% Adpated from [Halabi2009]
% Halabi, N., Rivoire, O., Leibler, S. & Ranganathan, R. 
% Protein sectors: evolutionary units of three-dimensional structure. 
% Cell 138, 774?86 (2009).


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

% Specifying the name of protein, and whether to run the shuffling code
protein = 'Nef';

%% Loading biochemical domains

biodomain = biochemical_domains(protein);

%% Preprocessing the data

load msa_nef
[B,Bcap,lambda,site_freq_no_mutation,true_indices,freq_bin,prev_aa,N,M,ls] = preprocessing_nef(msa);

%% SCA

B = double(~B); %definition of binary MSA is opposite in SCA with 0 for mutant and 1 for most prevalent amino acid at a site

freq_bg=[.073 .025 .050 .061 .042 .072 .023 .053 .064 .089...
    .023 .043 .052 .040 .052 .073 .056 .063 .013 .033];

% Background probabilities for the prevalent amino acids:
freq_bg_bin=freq_bg(prev_aa);

%SCA weighting matrix
W=log(freq_bin.*(1-freq_bg_bin)./(freq_bg_bin.*(1-freq_bin)));

%SCA matrix:
freq_pairs_bin=B'*B/N;
C_bin=freq_pairs_bin-freq_bin'*freq_bin;
C_sca=(W'*W).*abs(C_bin);
[eigvec_sca,lambda_sca] = eig_sort(C_sca);

%Shuffling the MSA to obtain number of significant eigenvectors
N_samples=100; %No. of randomized samples
lambda_rnd=zeros(N_samples,1);

for s=1:N_samples
    B_rnd=zeros(N,M);
    
    for pos=1:M
        perm_seq=randperm(N);
        B_rnd(:,pos)=B(perm_seq(:),pos);
    end
    
    freq_pairs_bin_rnd=B_rnd'*B_rnd/N;
    C_bin_rnd=freq_pairs_bin_rnd-freq_bin'*freq_bin;
    C_rnd = (W'*W).*abs(C_bin_rnd);
    
    [eigvect_unsorted,lambda_unsorted]=eig(C_rnd);
    lambda_rnd(s)=max(diag(lambda_unsorted));
end

alpha = sum(lambda_sca>max(lambda_rnd))

%The above procedure always results in alpha = 1 for SCA. However, we set
%the value of alpha equal to the one obtained using the Pearson
%correlation matrix based method (PCA and SPCA). Note that this value was
%chosen by inspection in Halabi et al. 2009.

alpha = 4

%% Forming SCA sectors

epsilon = 1.5/sqrt(M); %instead of 2/sqrt(M) as results for this are better

%sector with indices according to mutating sites (M = 451)
sec_eig_sca = cell(1,alpha);
%sector with indices according to Gag (ls = 500)
sec_eig_sca_true = cell(1,alpha);   
%sector with indices according to Gag and 100% conserved sites incorporated
sec_eig_sca_incl_cs = cell(1,alpha);

for kk = 2:alpha+1 
    sec_eig_sca{kk-1} = find(abs(eigvec_sca(:,kk))>epsilon);
    sec_eig_sca_true{kk-1} = true_indices(sec_eig_sca{kk-1});
    sec_eig_sca_incl_cs{kk-1} = combining_sectors_with_conserved_sites(sec_eig_sca_true{kk-1},site_freq_no_mutation,1,ls);
end

%% Statistical signficance of biochemical association of SCA sectors

%Table 2
sec_asso_sca = zeros(1,length(biodomain)); %sector associated to a particular biochemical domain
pvalue_asso_sca = zeros(1,length(biodomain)); %statistical significance of the sector-biochemical domain association

for kk = 1:length(biodomain)
    [sec_asso_sca(kk),pvalue_asso_sca(kk)] = ...
        compute_association(biodomain(kk).sites,sec_eig_sca_incl_cs,sec_eig_sca_true,ls,M);
end

fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('Significance of inferred %s sectors using SCA\n',protein)
fprintf('-----------------------------------------------------------------------------------\n')

for kk = 1:length(biodomain)
    fprintf('Sector %d is associated with %s (P = %.2e).\n',...
        sec_asso_sca(kk),biodomain(kk).name,pvalue_asso_sca(kk));
end
