function compare_size_sectors(sec_eig_spca_incl_cs,sec_eig_pca_incl_cs,alpha)

% Code to compare the size (number of sites) in RoCA and PCA sectors
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sec_eig_spca_incl_cs  a cell consisting of sectors formed using the
%                           RoCA method and supplemented by the conserved
%                           sites
%     sec_eig_pca_incl_cs   a cell consisting of sectors formed using the
%                           PCA method and supplemented by the conserved
%                           sites
%     alpha                 number of significant eigenvectors

%% Setting the font type and size for output figure

set(0,'DefaultAxesFontSize',22)
set(0,'DefaultTextFontSize',22)
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')

%% Computing the length of sectors formed using both methods

length_sec_RoCA = cellfun(@length,sec_eig_spca_incl_cs);
length_sec_PCA = cellfun(@length,sec_eig_pca_incl_cs);

length_sec_PCA_zeropadded = zeros(1,alpha);
for kk = 1:length(length_sec_PCA)
    length_sec_PCA_zeropadded(kk) = length_sec_PCA(kk);
end

%% Generating figure

figure;
bbb = bar(1:alpha,...
    [length_sec_PCA_zeropadded.' length_sec_RoCA.']);
legend('PCA','RoCA')
ylabel('Number of sites')
legend boxoff
xlim([0.5 alpha+.5])
ylim([0 200])
bbb(1).FaceColor = [228,26,28]/255;
bbb(2).FaceColor = [55,126,184]/255;
xlabel('Sector')

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'YTick'       , 0:40:200, ...
  'LineWidth'   , 1        );