function fig_boxplot(data,ITER,ylimdata,ytickdata,ylabeldata)

% Code to generate boxplot figure
% 
% Last updated: 2017-03-25

%%

whisker_value = 1.5;
widths_box_value = 0.5;

G = [ones(1,ITER) 2*ones(1,ITER) 3*ones(1,ITER) 4*ones(1,ITER)...
    5*ones(1,ITER) 6*ones(1,ITER) 7*ones(1,ITER) 8*ones(1,ITER)];

figure;
bh = boxplot(data,G,...
    'whisker',whisker_value,'symbol','ko',...
    'color','kkwkkwkk','medianstyle','target',...
    'labels',{'2','4','8','','','','',''},'widths',widths_box_value);
set(bh,'linewidth',2);
ylabel(ylabeldata)
ylim(ylimdata)
xlabel('N/M')
h=findobj(gca,'tag','Outliers');delete(h);%removing outliers for clean figures

%
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'YTick'       , ytickdata, ...
  'XTick'       , [1.5 4.5 7.5],...
  'LineWidth'   , 2        );
%

color_box{1} = [55,126,184]/255;
color_box{2} = [228,26,28]/255;
color_box{3} = 'w';
color_box{4} = [55,126,184]/255;
color_box{5} = [228,26,28]/255;
color_box{6} = 'w';
color_box{7} = [55,126,184]/255;
color_box{8} = [228,26,28]/255;
facealpha = [.6 .6 1 .6 .6 1 .6 .6];
h = findobj(gca,'Tag','Box');
len_excl_middle = setdiff(1:8,[3 6]);
for k=1:length(len_excl_middle)
    j = len_excl_middle(k);
    patch(get(h(j),'XData'),get(h(j),'YData'),color_box{j},'FaceAlpha',facealpha(j));
end
h=findobj(gca,'tag','MedianInner');
set(h,'Marker','o','MarkerSize',8,'MarkerFaceColor','none')
delete(h([3 6]))
h=findobj(gca,'tag','MedianOuter');
set(h,'Marker','o','MarkerSize',8,'MarkerFaceColor','none')
delete(h([3 6]))