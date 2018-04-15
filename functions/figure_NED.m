function figure_NED(alpha,NED_abs_random,NED_abs_inter,NED_abs_intra)

% Code to calculate Normalized Entropy Deviation (NED) for characterizing
% statistical independence between pairs of inferred sectors
% 
% Last updated: 2017-03-25
% 
% INPUT
%     alpha                 number of significant eigenvectors
%     NED_abs_random        NED for the randomized (null) case 
%     NED_abs_inter         Inter-sector NED 
%     NED_abs_intra         Intra-sector NED

%% Number of total pairs for inferred sectors

no_of_sec_pairs = nchoosek(alpha,2);

%% Increasing the smallest value of NED-Inter and NED-Random for display purposes

mean_NED_abs_random = mean(NED_abs_random);
small_value = 0.003;
for kk = 1:no_of_sec_pairs
    if mean_NED_abs_random(kk)<small_value
        mean_NED_abs_random(kk) = small_value;
    end
    if NED_abs_inter(kk)<small_value
        NED_abs_inter(kk) = small_value;
    end
end

%% Generating figure

figure('name','Statistical independence of SPCA sectors using Normalized Entropy Deviation (NED)',...
    'units','normalized','outerposition',[0 0 1 1])
pp = nchoosek(1:alpha,2); %possible pairs

data = [];
for kk = 1:no_of_sec_pairs
    data = [data; mean_NED_abs_random(kk) NED_abs_inter(kk) max(NED_abs_intra(pp(kk,1)),NED_abs_intra(pp(kk,2)))];
end
bbb = bar(1:no_of_sec_pairs,data);
bbb(1).FaceColor = [106,61,154]/255;%purple [228,26,28]/255;     %red
bbb(2).FaceColor = [150 150 150]./255;  %gray
bbb(3).FaceColor = [35,139,69]./255;%dark green, [55,126,184]./256;   %blue
set(gca, 'TickDir'     , 'out'     , ...
    'box', 'off',...
    'TickLength'  , [.01 .01])
legend('Random','Inter-sector','Intra-sector')
legend boxoff
axis([0.5 no_of_sec_pairs+.5 0 max(NED_abs_intra)])
set(gca,'XTick',1:no_of_sec_pairs)
if no_of_sec_pairs == 15 %Gag
    set(gca,'XTickLabel',{'(1,2)','(1,3)','(1,4)','(1,5)','(1,6)','(2,3)','(2,4)','(2,5)','(2,6)','(3,4)','(3,5)','(3,6)','(4,5)','(4,6)','(5,6)'})
else
    set(gca,'XTickLabel',{'(1,2)','(1,3)','(1,4)','(2,3)','(2,4)','(3,4)'})
end
% ylabel('Normalized Entropy Deviation (NED)')
ylabel('NED')
xlabel('Pairs of sectors')
