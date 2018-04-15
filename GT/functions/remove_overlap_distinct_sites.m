function [sec_eig,indicesSitesToBeDeleted,ties_unique] = ...
    remove_overlap_distinct_sites(sec_eig,distinct_sites,ties,C_clean)

%Code to assign overlapping site to the PCA sector with largest mean
%absolute correlation coefficient with respect to the distinct sites in 
%each sector
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sec_eig                   Initial sectors
%     distinct_sites            Distinct sites in each sector
%     ties                      Ties in each sector
%     C_clean                   Cleaned correlation matrix
%
% OUTPUT
%     sec_eig                   Output sectors with ties decided
%     indicesSitesToBeDeleted   Sites deleted from each sector
%     ties_unique               Unique ties

%%

sitesToBeDeleted = cell(1, length(sec_eig));
indicesSitesToBeDeleted = cell(1, length(sec_eig));

ties_unique = [];
for m = 1:length(sec_eig)
    ties_unique = [ties_unique; ties{m}];
end
ties_unique = unique(ties_unique); %Total ties

if sum(cellfun(@isempty,distinct_sites))~=0
    fprintf('\n Overlap cannot be removed as at least one sector has no distinct sites \n')
else
    for m = 1:length(ties_unique)
        mean_abs_coeff_for_this_site = zeros(1,length(sec_eig));
        for k = 1:length(sec_eig)
            mean_abs_coeff_for_this_site(k) = mean(abs(C_clean(ties_unique(m),distinct_sites{k})));
        end
        
        [~,indx_max_mean_abs_coeff_for_this_site] = ...
            max(mean_abs_coeff_for_this_site);
        
        removeSitefromSecs = setdiff(1:length(sec_eig),indx_max_mean_abs_coeff_for_this_site);
        
        for k = 1:length(removeSitefromSecs)
            sitesToBeDeleted{removeSitefromSecs(k)} = [sitesToBeDeleted{removeSitefromSecs(k)} ties_unique(m)];
            indicesSitesToBeDeleted{removeSitefromSecs(k)} = [indicesSitesToBeDeleted{removeSitefromSecs(k)} ...
                find(ismember(sec_eig{removeSitefromSecs(k)},ties_unique(m)))];
        end        
        
    end
    
    %Removing sites from sectors
    
    for k = 1:length(sec_eig)
        
        sec_eig{k}(indicesSitesToBeDeleted{k}) = [];
        
    end
    
end
