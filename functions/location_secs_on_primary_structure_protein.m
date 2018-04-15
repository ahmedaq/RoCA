function location_secs_on_primary_structure_protein(sec_eig_spca,sec_eig_pca,true_indices,freq_bin,M)

%Code superimposing both RoCA and RoCA sectors on primary structure of
%protein
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sec_eig_spca          a cell consisting of sectors formed using the
%                           RoCA method (indices according to M sites)
%     sec_eig_pca           a cell consisting of sectors formed using the
%                           RoCA method (indices according to M sites)
%     true_indices          correct indices of sites (according to the
%                           primary sequence)
%     freq_bin              1-by-M vector consisting of frequncy of
%                           dominant amino acid at each site
%     M                     total number of sites in MSA

%% Generating figure for PCA sectors

for kk = 1:length(sec_eig_pca)
    sec_pca(kk).def = sec_eig_pca{kk};
end

%color theme from qualitative set1-8
sec_pca(1).col = [228,26,28]/255;
sec_pca(2).col = [55,126,184]/255;
sec_pca(3).col = [77,175,74]/255;
sec_pca(4).col = [152,78,163]/255;
sec_pca(5).col = [255,127,0]/255;
sec_pca(6).col = [247,129,191]/255;

correct_primary_index = cell(1,M);
for kk = 1:M
    correct_primary_index{kk} = true_indices(kk);
end

figure('name','Location of sectors on primary structure of protein','units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
bar_secs_primary(-log(1-freq_bin),sec_pca,correct_primary_index)
axis([0 M+1 0 ceil(max(-log(1-freq_bin)))])
set(gca,'XTickLabel','');
set(gca, 'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01],...
  'box','off')
set(gca,'YTick',[0 ceil(max(-log(1-freq_bin)))]);
set(gca,'YTickLabel',[0 ceil(max(-log(1-freq_bin)))]);
title('PCA sectors','FontWeight','Normal')

%% Generating figure for RoCA sectors

for kk = 1:length(sec_eig_spca)
    sec_spca(kk).def = sec_eig_spca{kk};
end

%color theme from qualitative set1-8
sec_spca(1).col = [228,26,28]/255;
sec_spca(2).col = [55,126,184]/255;
sec_spca(3).col = [77,175,74]/255;
sec_spca(4).col = [152,78,163]/255;
sec_spca(5).col = [255,127,0]/255;
sec_spca(6).col = [247,129,191]/255;

correct_primary_index = cell(1,M);
for kk = 1:M
    correct_primary_index{kk} = true_indices(kk);
end

subplot(2,1,2)
bar_secs_primary(-log(1-freq_bin),sec_spca,correct_primary_index)
axis([0 M+1 0 ceil(max(-log(1-freq_bin)))])
set(gca, 'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01],...
  'box','off')
set(gca,'YTick',[0 ceil(max(-log(1-freq_bin)))]);
set(gca,'YTickLabel',[0 ceil(max(-log(1-freq_bin)))]);
ylabel('$$-\log(f_i)$$','interpreter','latex')
xlabel('$$i$$','interpreter','latex')
title('RoCA sectors','FontWeight','Normal')
