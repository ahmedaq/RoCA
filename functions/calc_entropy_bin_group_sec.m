function Entropy_bin_group_sec = calc_entropy_bin_group_sec(sec,msa_bin)

% Code to calculate overall entropy of sector for all possible 
% combinations of sites in it
%
% Last updated: 2017-03-25
% 
% INPUT
%     sec                   Sector sites
%     msa_bin               Binarized MSA
%
% OUTPUT
%     Entropy_bin_group_sec Overall entropy

%%

Entropy_bin_group_sec = 0;
total_combs = get_binary_vecs(length(sec));
freq_comb = zeros(length(total_combs),1);

for kk = 1:size(total_combs,1)
    
    comb_curent = repmat(total_combs(kk,:),size(msa_bin,1),1);
    freq_comb(kk) = sum(sum(msa_bin(:,sec)==comb_curent,2)==length(sec))/size(msa_bin,1);
    
    if freq_comb(kk)>0
        Entropy_bin_group_sec = Entropy_bin_group_sec - freq_comb(kk).*log(freq_comb(kk));
        
    end
end
