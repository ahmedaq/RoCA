function pvalue = pvalue_fisher(total_no_of_sites,no_of_sites_in_sector,sites_test,sites_test_in_sector)

% Code to compute P value using Fisher's exact test
% 
% Last updated: 2017-03-25
% 
% INPUT
%     total_no_of_sites         Total number of sites in protein
%     no_of_sites_in_sector     Total number of sites in inferred sector
%     sites_test                Total number of biologically significant 
%                               sites
%     sites_test_in_sector      Number of biologically significant sites
%                               present in the inferred sector
%
% OUTPUT
%     pvalue                    Associated P value

%%
warning off

pvalue = zeros(1,length(sites_test));

for j = 1:length(sites_test)    
    for k = sites_test_in_sector(j):1:min(sites_test(j),no_of_sites_in_sector)
        pvalue(j) = pvalue(j) + ...
            (nchoosek(sites_test(j),k)*...
            nchoosek(total_no_of_sites-sites_test(j),no_of_sites_in_sector-k))/...
            nchoosek(total_no_of_sites,no_of_sites_in_sector);
    end
end