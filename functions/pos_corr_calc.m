function [no_of_pairs_with_pos_corr,per_pos_corr,sum_pos_corr,sites_pos_corr_corr]...
    = pos_corr_calc(sector,corr_matrix,pos_thresh)

% Code to obtain percentage of positively correlated sites in a given
% sector
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sector                    The sector consisting of sites whose 
%                               positive correlation is to be calculated
%     corr_matrix               The correlation matrix whose entries are 
%                               to be checked
%     pos_thresh                The positive correlation threshold - The 
%                               sites are checked for correlations greater
%                               than this value
%
% OUTPUT
%     no_of_pairs_with_pos_corr Number of sites with correlation greater 
%                               than pos_thresh
%     per_pos_corr              Percentage of positively correlaled pairs 
%                               in this sector
%     sum_pos_corr              Sum of positive correlations of sites with
%                               corr>pos_thresh
%     per_pos_corr_sites        Percentage of 'unique' positively 
%                               correlated sites in this sector

%%

if nargout <= 3
    
    no_of_pairs_with_pos_corr = 0;
    sum_pos_corr = 0;
    for i = 1:length(sector)
        for j = i:length(sector)
            m = sector(i);
            n = sector(j);
            if corr_matrix(m,n) >= pos_thresh && m~=n
                no_of_pairs_with_pos_corr = no_of_pairs_with_pos_corr+1;
                sum_pos_corr = sum_pos_corr + corr_matrix(m,n);
            end
        end
    end
    
    per_pos_corr = no_of_pairs_with_pos_corr/(length(sector)*(length(sector)-1)/2)*100;
    
else
    
    no_of_pairs_with_pos_corr = 0;
    sum_pos_corr = 0;
    sites_pos_corr_corr = [];
    for i = 1:length(sector)
        for j = i:length(sector)
            m = sector(i);
            n = sector(j);
            if corr_matrix(m,n) >= pos_thresh && m~=n
                no_of_pairs_with_pos_corr = no_of_pairs_with_pos_corr+1;
                sum_pos_corr = sum_pos_corr + corr_matrix(m,n);
                sites_pos_corr_corr = [sites_pos_corr_corr; m n corr_matrix(m,n)];
            end
        end
    end
    
    per_pos_corr = no_of_pairs_with_pos_corr/(length(sector)*(length(sector)-1)/2)*100;
end