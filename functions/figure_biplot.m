function figure_biplot(sector,eigvec1,eigvec2,color_theme,M,n_secs,markersize,jitter)

% Code to generate each subplot in the biplot figure
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sector                A cell consisting of sectors formed using the
%                           PCA/RoCA method (indices according to M sites)
%     eigvec1               First eigenvector
%     eigvec2               Second eigenvector
%     color_theme           Color theme for sites in different sectors
%     M                     total number of sites in MSA
%     n_secs                Number of sectors formed using PCA/RoCA method
%     markersize            The size of circles in the biplots
%     jitter                The standard deviation of random noise added to
%                           each point in the biplot to avoid super-
%                           imposition of data points

%% Finding the overlapping and non-overlapping sites in sectors
sec_all=[];
for kkk = 1:n_secs
    sec_all = [sec_all;sector{kkk}];
end

rest = setdiff(1:M,unique(sec_all));

sites_overlap = [];
for kkk = 1:n_secs
    for mmm = kkk+1:n_secs
        sites_overlap = [sites_overlap ; intersect(sector{kkk},sector{mmm})];
    end
end

for kkk = 1:n_secs
    sector{kkk} = setdiff(sector{kkk},sites_overlap); %only non overlapping sites.
end

%% Generating figure

if isempty(sites_overlap)==1
    
    for kkk = 1:n_secs
        plot(eigvec1(sector{kkk}),eigvec2(sector{kkk}),'o','MarkerFaceColor',color_theme(kkk,:),...
            'MarkerEdgeColor','k','MarkerSize',markersize)
        hold on
    end
    plot(eigvec1(rest),eigvec2(rest),'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',markersize)
else
    
    for kkk = 1:n_secs
        plot(eigvec1(sector{kkk})+jitter*randn(length(sector{kkk}),1),eigvec2(sector{kkk})+jitter*randn(length(sector{kkk}),1),...
            'o','MarkerFaceColor',color_theme(kkk,:),'MarkerEdgeColor','k','MarkerSize',markersize)
        hold on
    end
    grey_color = 190*ones(1,3)/255;
    plot(eigvec1(sites_overlap),eigvec2(sites_overlap),'o',...
        'MarkerFaceColor',grey_color,'MarkerEdgeColor','k','MarkerSize',markersize-2)    
    plot(eigvec1(rest),eigvec2(rest),'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',markersize)
      
end

axis([round2dp(min(eigvec1),2)-.05 round2dp(max(eigvec1),2)+.05 ...
    round2dp(min(eigvec2),2)-.05 round2dp(max(eigvec2),2)+.05])

set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'XColor'      , [0 0 0], ...
    'YColor'      , [0 0 0], ...
    'YAxisLocation', 'left',...
    'XAxisLocation', 'top',...
    'LineWidth'   , 1    );

% 'YTick'       , round2dp(min(eigvec2),1):.1:round2dp(max(eigvec2),1),...
%     'XTick'       , round2dp(min(eigvec1),1):.1:round2dp(max(eigvec1),1),...