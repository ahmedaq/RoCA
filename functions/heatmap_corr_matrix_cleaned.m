function per_sites_in_sectors = heatmap_corr_matrix_cleaned(C_hat,sec_eig,sec_eig_incl_cs,alpha,ls,method)

%Code to generate heatmap of the cleaned correlation matrices (restricted
%to only the sector sites or complete one)
% 
% Last updated: 2017-03-25
% 
% INPUT
%     C_hat                 M-by-M cleaned correlation matrix 
%     sec_eig               a cell consisting of sectors formed using the
%                           PCA method (indices according to M sites)
%     sec_eig_incl_cs       a cell consisting of sectors formed using the
%                           PCA method and supplemented by the conserved
%                           sites
%     alpha                 number of significant eigenvectors
%     ls                    actual number of sites in the protein (before
%                           preprocessing)
%     method                options: 'RoCA' or 'PCA'
%
% OUTPUT
%     per_sites_in_sectors  total percentage of sites in sectors (including
%                           the conserved sites) respect to actual protein
%                           length (ls)


%% If method is not specified, generate results for RoCA

if nargin < 6
    method = 'RoCA';
end

%% Ordering the sites in sectors
M = length(C_hat);
sec_eig_ord = cell(1,alpha);

for kk = 1:alpha
   sec_eig_ord{kk} = order_sectors_according_C(sec_eig{kk},C_hat);
end

%% Arranging the sites according to sectors in a vector 'secs_cascaded_ord'
%% and assigning an overlapping site to only one sector

secs_cascaded_ord = [];
for kk = 1:alpha
    %Removing any overlapping site in subsequent sectors
    if sum(ismember(secs_cascaded_ord,sec_eig_ord{kk}))~=0
        sec_eig_ord{kk}((ismember(sec_eig_ord{kk},secs_cascaded_ord)))=[];
    end
    secs_cascaded_ord = [secs_cascaded_ord;sec_eig_ord{kk}];
end

%% Percentage of sites in all sectors (including the conserved sites)
secs_cascaded_incl_cs_ord = [];
sec_eig_incl_cs2 = sec_eig_incl_cs;
for m = 1:alpha
    %Removing any overlapping site in subsequent sectors
    if sum(ismember(secs_cascaded_incl_cs_ord,sec_eig_incl_cs2{m}))~=0
        sec_eig_incl_cs2{m}((ismember(sec_eig_incl_cs2{m},secs_cascaded_incl_cs_ord)))=[];
    end
    secs_cascaded_incl_cs_ord = [secs_cascaded_incl_cs_ord sec_eig_incl_cs2{m}];
    
end

per_sites_in_sectors = length(secs_cascaded_incl_cs_ord)/ls*100;

%% Generating heat map figures

if strcmp(method,'RoCA')
    % Heat map restricted to sites in sectors (Figure 4)
    heatmap_corr_matrix(C_hat(secs_cascaded_ord,secs_cascaded_ord),'northoutside',...
        sprintf('%s - Cleaned correlation matrix (only sectors)',method))
    
    %Complete heat map (Figure 1c)    
    heatmap_corr_matrix(C_hat([secs_cascaded_ord.' setdiff(1:M,unique(secs_cascaded_ord))],...
        [secs_cascaded_ord.' setdiff(1:M,unique(secs_cascaded_ord))]),'southoutside',...
        sprintf('%s - Cleaned correlation matrix (complete)',method))
    
else
    
    %Complete heat map (Figure 1c)
    heatmap_corr_matrix(C_hat([secs_cascaded_ord.' setdiff(1:M,unique(secs_cascaded_ord))],...
        [secs_cascaded_ord.' setdiff(1:M,unique(secs_cascaded_ord))]),'southoutside',...
        sprintf('%s - Cleaned correlation matrix (complete)',method))
    
end