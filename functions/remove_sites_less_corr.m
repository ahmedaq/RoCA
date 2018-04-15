function [sector_out] = ...
    remove_sites_less_corr(sector_in,sites_in_all_sectors,C_clean,pos_thresh,neg_thresh)

%Code to remove sites from a PCA sector which are not significantly
%correlated with the remaining sites in the sector
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sector_in             one PCA sector
%     sites_in_all_sectors  same PCA sector or sites in all sectors
%     C_clean               cleaned correlation matrix
%     pos_thresh            threshold defining significnat positive
%                           correlations
%     neg_thresh            threshold defining significnat negative
%                           correlations
%
% OUTPUT
%     sector_out            sector with sites removed

%%

if length(sector_in)>1 %size of sector should be greater than at least 1
    
    sites_drop = [];
    indx_drop=[];
    for k = 1:length(sector_in)
        site_check = sector_in(k);
        temp = sites_in_all_sectors;
        temp((ismember(site_check,sites_in_all_sectors))) = [];  %Removing the site_check to avoid autocorr
        
        if neg_thresh ~= 0
            if max((C_clean(site_check,temp)))...
                    < pos_thresh && min((C_clean(site_check,temp))) > neg_thresh
                sites_drop = [sites_drop site_check];
                indx_drop = [indx_drop k];
            end
        else
            if max((C_clean(site_check,temp)))< pos_thresh
                sites_drop = [sites_drop site_check];
                indx_drop = [indx_drop k];
            end
        end
    end
    
    sector_out = sector_in;
    sector_out(indx_drop) = [];
    
else
    sector_out = sector_in;
    sites_drop=[];
    indx_drop = [];
end

