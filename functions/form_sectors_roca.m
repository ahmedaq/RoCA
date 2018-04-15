function [sec_eig, sec_eig_true, sec_eig_incl_cs, length_sec] = ...
    form_sectors_roca(PC,alpha,site_freq_no_mutation,true_indices,ls)

% Code to form sectors using the RoCA method
% 
% Last updated: 2017-03-25
% 
% INPUT
%     PC                    M-by-(alpha) matrix, columns corresponding to sparse 
%                           estimate of leading eigenvectors
%     alpha                 number of significant eigenvectors
%     site_freq_no_mutation indices of 100% conserved sites
%     true_indices          correct indices of sites (according to the
%                           primary sequence)
%     ls                    actual number of sites in the protein (before
%                           preprocessing)
%     
% OUTPUT
%     sec_eig               1-by-(alpha) cell consisting of sectors formed
%                           using the RoCA method (indices according to
%                           M sites)
%     sec_eig_true          1-by-(alpha) cell consisting of sectors formed
%                           using the RoCA method (indices according to 
%                           the correct ls sites)
%     sec_eig_incl_cs       1-by-(alpha) cell consisting of sectors formed
%                           using the RoCA method and supplemented by the 
%                           100% conserved sites
%     sec_eig               1-by-(alpha) vector consisting of length of 
%                           each RoCA sector

%%
M = length(PC(:,1));
%sector with indices according to mutating sites (M = 451)
sec_eig = cell(1,alpha);
%sector with indices according to Gag (ls = 500)
sec_eig_true = cell(1,alpha);   
%sector with indices according to Gag and 100% conserved sites incorporated
sec_eig_incl_cs = cell(1,alpha);
%length of each sector with 100% conserved sites incorporated
length_sec = zeros(1,alpha);    

for kk = 1:alpha
    sec_eig{kk} = find(abs(PC(:,kk))>1/sqrt(M));
    sec_eig_true{kk} = true_indices(sec_eig{kk});
    sec_eig_incl_cs{kk} = combining_sectors_with_conserved_sites(sec_eig_true{kk},site_freq_no_mutation,1,ls);
    length_sec(kk) = length(sec_eig_incl_cs{kk});
end
