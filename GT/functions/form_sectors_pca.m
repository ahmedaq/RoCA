function [Q_sample, sec_eig, sec_eig_true, sec_eig_incl_cs, length_sec] = ...
    form_sectors_pca(msa,alpha,site_freq_no_mutation,true_indices,ls,pos_thresh,neg_thresh)

% Code to form sectors using the PCA method proposed in Quadeer et al. 2014
% 
% Last updated: 2017-03-25
% 
% INPUT
%     msa                   N-by-M MSA matrix
%     alpha                 number of significant eigenvectors
%     site_freq_no_mutation indices of 100% conserved sites
%     true_indices          correct indices of sites (according to the
%                           primary sequence)
%     ls                    actual number of sites in the protein (before
%                           preprocessing)
%     pos_thresh            threshold defining significnat positive
%                           correlations
%     neg_thresh            threshold defining significnat negative
%                           correlations
%     
% OUTPUT
%     Q_sample              M-by-(alpha) matrix, columns corresponding to
%                           the significant eigenvectors
%     sec_eig               a cell consisting of sectors formed using the
%                           PCA method (indices according to M sites)
%     sec_eig_true          a cell consisting of sectors formed using the
%                           PCA method (indices according to the correct 
%                           ls sites)
%     sec_eig_incl_cs       a cell consisting of sectors formed using the
%                           PCA method and supplemented by the conserved
%                           sites
%     length_sec            a vector consisting of length of each PCA
%                           sector


[~,M]=size(msa);

%% Obtaining sample eigenvectors and corresponding cleaned correlation matrix

C = corrcoef(msa);
[Q_sample,lambda_sample] = eig_sort(C);

C_clean = zeros(M,M);
for kk = 1:alpha
    C_clean = C_clean + Q_sample(:,kk)*lambda_sample(kk)*(Q_sample(:,kk)');
end


%% Initial sectors
epsilon = 0.005;
for kk = 1:alpha
    sec_eig{kk} = find(abs(Q_sample(:,kk))>epsilon);
end

%% Postprocessing I: 
%% Removing overlap in sectors

for i = 1:alpha
    secs_rest = setdiff(1:alpha,i);
    sites_rest = [];
    for j = 1:length(secs_rest)
        sites_rest = [sites_rest sec_eig{secs_rest(j)}.'];
    end
    distinct_sites{i} = setdiff(sec_eig{i},sites_rest);
    ties{i} = sec_eig{i}(ismember(sec_eig{i},sites_rest));
end

if min(cellfun(@length,distinct_sites))>1 %more than one distinct site in a sector
    [sec_eig] = remove_overlap_distinct_sites(sec_eig,distinct_sites,ties,C_clean);
    
else
    for m = 1:alpha
        for k = m+1:alpha
            [sec_eig{m},sec_eig{k}] = remove_overlap_sectors(sec_eig{m},sec_eig{k},C_clean);
        end
    end
    
end

%% Postprocessing II: 
%% Removing sites with insignificant correlations (between negative and positive) threshold from the sectors

for m = 1:alpha
    sec_eig{m} = ...
        remove_sites_less_corr(sec_eig{m},sec_eig{m},C_clean,pos_thresh,neg_thresh);
end

%% Removing empty sectors

sec_eig(cellfun(@isempty,sec_eig))=[];

%% Forming sectors with true indices and incorporating 100% conserved sites in the sectors

%sector with indices according to Gag (ls = 500)
sec_eig_true = cell(1,length(sec_eig));   
%sector with indices according to Gag and 100% conserved sites incorporated
sec_eig_incl_cs = cell(1,length(sec_eig));
%length of each sector with 100% conserved sites incorporated
length_sec = zeros(1,length(sec_eig));   

for kk = 1:length(sec_eig)
    sec_eig_true{kk} = true_indices(sec_eig{kk});
    sec_eig_incl_cs{kk} = combining_sectors_with_conserved_sites(sec_eig_true{kk},site_freq_no_mutation,1,ls);
    length_sec(kk) = length(sec_eig_incl_cs{kk});
end