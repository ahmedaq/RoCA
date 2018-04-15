function [no_of_pairs_with_neg_corr,per_neg_corr,sum_neg_corr,per_neg_corr_sites,sites_neg_corr_corr] = ...
    neg_corr_calc(sector,corr_matrix,neg_thresh)

% Code to obtain percentage of negatively correlated sites in a given
% sector
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sector                    The sector consisting of sites whose 
%                               negative correlation is to be calculated
%     corr_matrix               The correlation matrix whose entries are 
%                               to be checked
%     neg_thresh                The negative correlation threshold - The 
%                               sites are checked for correlations smaller
%                               than this value
%
% OUTPUT
%     no_of_pairs_with_neg_corr Number of sites with correlation smaller 
%                               than neg_thresh
%     per_neg_corr              Percentage of negatively correlaled pairs 
%                               in this sector
%     sum_neg_corr              Sum of negative correlations of sites with
%                               corr<neg_thresh
%     per_neg_corr_sites        Percentage of 'unique' negatively 
%                               correlated sites in this sector

%%

if nargout <= 3
    
    no_of_pairs_with_neg_corr = 0;
    sum_neg_corr = 0;
    for i = 1:length(sector)
        for j = i:length(sector)
            m = sector(i);
            n = sector(j);
            if corr_matrix(m,n) < neg_thresh
                no_of_pairs_with_neg_corr = no_of_pairs_with_neg_corr+1;
                sum_neg_corr = sum_neg_corr + corr_matrix(m,n);
            end
        end
    end
    
    per_neg_corr = no_of_pairs_with_neg_corr/(length(sector)*(length(sector)-1)/2)*100;
    
else
    
    no_of_pairs_with_neg_corr = 0;
    sum_neg_corr = 0;
    sites_neg_corr_corr = [];
    for i = 1:length(sector)
        for j = i:length(sector)
            m = sector(i);
            n = sector(j);
            if corr_matrix(m,n) < neg_thresh
                no_of_pairs_with_neg_corr = no_of_pairs_with_neg_corr+1;
                sum_neg_corr = sum_neg_corr + corr_matrix(m,n);
                sites_neg_corr_corr = [sites_neg_corr_corr; m n corr_matrix(m,n)];
            end
        end
    end
    
    per_neg_corr = no_of_pairs_with_neg_corr/(length(sector)*(length(sector)-1)/2)*100;
    per_neg_corr_sites = length(unique(sites_neg_corr))/length(sector)*100;
end