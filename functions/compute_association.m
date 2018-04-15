function [sec_association,pvalue] = compute_association(bio_domain,sec_eig_incl_cs,sec_eig_true,ls,M)

%Code to compute the most significnat association of an experimentally known
%biochemical domain of the protein to an inferred sector
% 
% Last updated: 2017-03-25
% 
% INPUT
%     bio_domain            sites important for a particular biochemical
%                           domain
%     sec_eig_incl_cs       a cell consisting of sectors formed using the
%                           PCA/RoCA method and supplemented by the conserved
%                           sites
%     sec_eig_true          1-by-(alpha) cell consisting of sectors formed
%                           using the RoCA method (indices according to 
%                           the correct ls sites)
%     ls                    actual number of sites in the protein (before
%                           preprocessing)
%     M                     total number of sites in MSA
%
% OUTPUT
%     sec_association       Number of the sector associated with the input
%                           bio_domain
%     pvalue                Statistical significance (P-value) associated
%                           with the most significant association

%% Calculating the pvalue of association of all sectors to the input biodomain

for kk = 1:length(sec_eig_incl_cs)
    
    pvalue_association_incl_cs(kk) = ...
        pvalue_fisher(ls,length(sec_eig_incl_cs{kk}),...
        length(bio_domain),...
        sum(ismember(bio_domain,sec_eig_incl_cs{kk})));
    
    pvalue_association(kk) = ...
        pvalue_fisher(M,length(sec_eig_true{kk}),...
        length(bio_domain),...
        sum(ismember(bio_domain,sec_eig_true{kk})));
    
end

%% Finding the most significant associated sector

sec_association = find(pvalue_association_incl_cs<=0.05);

if length(sec_association)>1 %if multiple sectors are associated with P<0.05
    sec_association = sec_association(pvalue_association(sec_association)==min(pvalue_association(sec_association)));
elseif isempty(sec_association) %no significant association
    sec_association = find(pvalue_association_incl_cs==min(pvalue_association_incl_cs));
    sec_association = sec_association(1);
end
pvalue = pvalue_association_incl_cs(sec_association);
