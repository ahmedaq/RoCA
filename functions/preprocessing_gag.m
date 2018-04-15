function [B,Bcap,lambda,site_freq_no_mutation,true_indices,freq_bin,prev_aa,N,M,ls] = preprocessing_gag(A)

%Code for preprocessing the HIV Gag data
% 
% Last updated: 2017-03-25
% 
% INPUT
%     A                     Downloaded amino acid MSA 
%
% OUTPUT
%     B                     N-by-M matrix of binary coded MSA
%     Bcap                  N-by-M MSA matrix with effect of phylogeny removed
%     lambda                Sample eigenvalues of Bcap
%     site_freq_no_mutation Indices of conserved sites
%     true_indices          correct indices of sites (according to the
%                           primary sequence)
%     freq_bin              1-by-M vector consisting of frequency of
%                           dominant amino acid at each site          
%     prev_aa               1-by-M vector consisting of identity of 
%                           dominant amino acid at each site          
%     N                     Number of sequences after preprocessing
%     M                     Number of sites after preprocessing
%     ls                    actual number of sites in the protein (before
%                           preprocessing)

%% Identifying and removing clusters of outlying sequences

outlying_seqs = obtain_outlying_seqs_gag(A);
A(outlying_seqs,:) = [];

%% Removing the 100% conserved sites

[ns,ls]=size(A); %ns=no.of seq, ls=length of seq
profile_seq = seqprofile(A);
[~,pos_wt] = max(profile_seq);
Cseq = int2aa(pos_wt);
Cseq_mtrx = repmat(Cseq,ns,1);
msa_binary = double(A~=Cseq_mtrx);
%replacing mutation by 1 and no mutation by 0

freq_single_mutation = mean(msa_binary); %frequency of mutation at each site 
site_freq_no_mutation = find(freq_single_mutation==0); %100% conserved sites

%% Finding sites with only one mutation and gap/blank/ambiguous amino acid on it

no_of_mutations_per_site = sum(msa_binary,1); %no of mutations at each site
sites_with_one_mutation = find((no_of_mutations_per_site)==1);
sites_with_only_blank_mutation = [];
mutants_sites_with_one_mutation = cell(1,length(sites_with_one_mutation));
length_mutants_sites_with_one_mutation = zeros(1,length(sites_with_one_mutation));
for kk = 1:length(sites_with_one_mutation)
    mutants_sites_with_one_mutation{kk} = A(A(:,sites_with_one_mutation(kk))~=Cseq(sites_with_one_mutation(kk)),sites_with_one_mutation(kk));
    length_mutants_sites_with_one_mutation(kk) = length(mutants_sites_with_one_mutation{kk});
    if mutants_sites_with_one_mutation{kk} == 'B'
        sites_with_only_blank_mutation= [sites_with_only_blank_mutation sites_with_one_mutation(kk)];
    end
end

%%  Finding sites with more than 12% gaps/blanks/ambiguous amino acids
no_of_blanks = zeros(1,ls);
for m = 1:ls
    temp = 0;
    for k = 1:ns
        if A(k,m) == 'B'
            temp = temp +1;
        end
    end
    no_of_blanks(m) = temp;
end
per_no_of_blanks = no_of_blanks/ns*100;
site_blanks_gt_12pc = find(per_no_of_blanks>12);

%% Combining all sites to be removed
site_freq_no_mutation = unique([site_freq_no_mutation sites_with_only_blank_mutation site_blanks_gt_12pc]);
percentage_conserved_sites = length(site_freq_no_mutation)/ls*100;

%variable to save primary structure indices of the sites remaining in MSA
true_indices = setdiff([1:ls],site_freq_no_mutation);

%% Removing all the conserved and problematic sites

A(:,site_freq_no_mutation)=[];

fprintf('Data statistics')
fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('Percentage of conserved/problematic sites = %.2f.\n',percentage_conserved_sites);

[N,M] = size(A);

fprintf('Number of sequences after data preprocessing, N = %d,\nNumber of sites after data preprocessing, M = %d.\n',N,M);

%% Binary approximation and its validation

[B,freq_bin,prev_aa] = binary_approx(A);

%% Removing phylogeny

[Bcap,lambda] = remove_phylogeny(B,1);