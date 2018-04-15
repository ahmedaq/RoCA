function [sec1_new,sec2_new,sites_del_sec1,sites_del_sec2] = ...
    remove_overlap_sectors(sec1,sec2,C)

%Code to assign overlapping site to the PCA sector with largest mean
%absolute correlation coefficient with the sites in that sector
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sec1                  sector 1
%     sec2                  sector 2
%     C                     Correlation matrix
%
% OUTPUT
%     sec1_new              sector 1 after deciding overlapping sites
%     sec2_new              sector 2 after deciding overlapping sites
%     sites_del_sec1        sites removed from sector 1
%     sites_del_sec2        sites removed from sector 2

%%

indx_site_S1_present_in_S2 = find(ismember(sec1,sec2));
site_S1_present_in_S2 = sec1(indx_site_S1_present_in_S2);
sites_del_S1_comp_S2 = [];
indx_del_S1_comp_S2 = [];
sites_del_S2_comp_S1 = [];

for k = 1:length(site_S1_present_in_S2)
    if mean(abs(C(site_S1_present_in_S2(k),sec1))) < ...
            mean(abs(C(site_S1_present_in_S2(k),sec2)))
        sites_del_S1_comp_S2 = [sites_del_S1_comp_S2 site_S1_present_in_S2(k)];
        indx_del_S1_comp_S2 = [indx_del_S1_comp_S2 indx_site_S1_present_in_S2(k)];
    else
        sites_del_S2_comp_S1 = [sites_del_S2_comp_S1 site_S1_present_in_S2(k)];        
    end    
end

sec1_new = sec1;
sec1_new(indx_del_S1_comp_S2)=[];    %Removing sites
sec2_new = sec2;
sec2_new(ismember(sec2_new,sites_del_S2_comp_S1))=[]; %Removing sites
sites_del_sec1 = sites_del_S1_comp_S2;
sites_del_sec2 = sites_del_S2_comp_S1;