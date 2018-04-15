function generate_biplots(PC,sec_eig,alpha,markersize,jitter,method)

% Code to generate the biplots of the eigenvectors
% 
% Last updated: 2017-03-25
% 
% INPUT
%     PC                    M-by-(alpha) matrix, columns corresponding to sparse 
%                           estimate of leading eigenvectors
%     sec_eig               a cell consisting of sectors formed using the
%                           PCA/RoCA method (indices according to M sites)
%     alpha                 number of significant eigenvectors
%     markersize            the size of circles in the biplots
%     jitter                the standard deviation of random noise added to
%                           each point in the biplot to avoid super-
%                           imposition of data points
%     method                options: 'RoCA' or 'PCA'

%% If method is not specified, generate results for RoCA

if nargin < 6
    method = 'RoCA';
end

%% Setting the color theme for biplots

%color theme from qualitative set1-8
color_theme(1,:)=[228,26,28]/255;
color_theme(2,:)=[55,126,184]/255;
color_theme(3,:)=[77,175,74]/255;
color_theme(4,:)=[152,78,163]/255;
color_theme(5,:)=[255,127,0]/255;
color_theme(6,:)=[247,129,191]/255;

%% Generating figure

M = length(PC(:,1));

figure('name',sprintf('Biplots of %s eigenvectors',method),'units','normalized','outerposition',[0 0 1 1])
for kk = 1:alpha-1
    count = (kk-1)*alpha+1;
    
    for mm = kk+1:alpha
        subplot(alpha-1,alpha-1,count)
        
        figure_biplot(sec_eig,PC(:,mm),PC(:,kk),...
            color_theme,M,length(sec_eig),markersize,jitter)
        
        count = count+1;
        
        if kk == 1
            xlabel(sprintf('PC %d',mm),'FontSize',22)
        end
        
        if mm == kk+1
            ylabel(sprintf('PC %d',kk),'FontSize',22)
        end
        
    end
end


