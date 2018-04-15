function overlap = calculate_overlap(sec_eig_spca,sec_eig_pca,n_secs_spca,n_secs_pca)

% Code to calculate and plot the merging of RoCA sectors in the PCA ones
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sec_eig_spca          a cell consisting of sectors formed using the
%                           RoCA method (indices according to M sites)
%     sec_eig_pca           a cell consisting of sectors formed using the
%                           RoCA method (indices according to M sites)
%     n_secs_spca           Number of sectors formed using the RoCA method
%     n_secs_pca            Number of sectors formed using the PCA method

%% Defining color for the bar plot

color_darkgrey = [166 86 40]/255;

%% Calculating overlap

for kk = 1:n_secs_pca
    for mm = 1:n_secs_spca
        overlap(kk,mm) = sum(ismember(sec_eig_spca{mm},sec_eig_pca{kk}))/length(sec_eig_spca{mm})*100;
    end
end

%% Generating figure

figure('name','Percentage overlap of sectors','units','normalized','outerposition',[0 0 1 1])
for kk = 1:n_secs_pca
    subplot(ceil(n_secs_pca/2),2,kk)
    bar(1:n_secs_spca,overlap(kk,:),0.3,'FaceColor',color_darkgrey,'EdgeColor','k')
    ylabel('Percentage overlap')
    xlabel('RoCA sector')
    title(sprintf('PCA sector %d',kk))
    axis([0.5 n_secs_spca+0.5 0 100])
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'on'      , ...
        'XGrid'       , 'off'      , ...
        'YGrid'       , 'on'      , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'YTick'       , 0:25:100, ...
        'LineWidth'   , 1        );
end

