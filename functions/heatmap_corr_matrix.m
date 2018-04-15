function heatmap_corr_matrix(C,colorbar_position,title)

%Code to generate heatmap of the input correlation matrix
% 
% Last updated: 2017-03-25
% 
% INPUT
%     C                     Input correlation (square) matrix
%     colorbar_position     Position of the colorbar in the figure
%     title                 Title of the figure

%% Generating figure
warning off
figure('name',title)
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1); 

% map = [103,0,31;
% 214,96,77;
% 253,219,199;
% 247,247,247;
% 209,229,240;
% 146,197,222;
% 67,147,195;
% 33,102,172;
% 5,48,97]./255;

% map = brewermap(200,'RdBu');
% map2 = [map(1:2:100,:); map(101:end,:)];

imshow((C),'Colormap',jet,'DisplayRange',[-0.2 .4]);

h = colorbar(colorbar_position);

set(h,...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , .02 , ...
    'Ytick',[-.2 0 .2 .4])

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [0 0] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , 0:0.1:0.1, ...
  'LineWidth'   , 1        );

