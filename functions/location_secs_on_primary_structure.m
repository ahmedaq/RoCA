function location_secs_on_primary_structure(sec_eig,true_indices,freq_bin,M)

%Code to generate bar plot of sector sites on primary structure
% 
% INPUT
%     sec_eig               1-by-(alpha) cell consisting of sectors formed
%                           using the SPCA method (indices according to
%                           M sites)
%     true_indices          correct indices of sites (according to the
%                           primary sequence)    
%     freq_bin              1-by-M vector consisting of frequency of
%                           dominant amino acid at each site
%     M                     total number of sites in MSA

%% Defining sector sites and color

for kk = 1:length(sec_eig)
    sec(kk).def = sec_eig{kk};
end

%color theme from qualitative set1-8
sec(1).col = [228,26,28]/255;
sec(2).col = [55,126,184]/255;
sec(3).col = [77,175,74]/255;
sec(4).col = [152,78,163]/255;
sec(5).col = [255,127,0]/255;
sec(6).col = [247,129,191]/255;

%% Making a cell of correct indices for each site

for kk = 1:M
    correct_primary_index{kk} = true_indices(kk);
end

%% Generating figure

figure('name','Location of sectors on primary structure of protein','units','normalized','outerposition',[0 0 1 1])
bar_secs_primary(-log(1-freq_bin),sec,correct_primary_index)
axis([0 M+1 0 ceil(max(-log(1-freq_bin)))])
set(gca, 'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01],...
  'box','off')
set(gca,'YTick',[0 ceil(max(-log(1-freq_bin)))]);
set(gca,'YTickLabel',[0 ceil(max(-log(1-freq_bin)))]);
ylabel('$$-\log(f_i)$$','interpreter','latex')
xlabel('$$i$$','interpreter','latex')
title('Sectors','FontWeight','Normal')

