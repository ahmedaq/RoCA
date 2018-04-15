function bar_secs_primary(V,sec,correct_index)

%Code to generate bar plot of sector sites on primary structure
% 
% (adapted from Halabi et al. [1])
% [1] Halabi, N. et al. Supplemental Data Theory Protein Sectors: 
%     Evolutionary Units of Three-Dimensional Structure. Cell 138, (2009).
% 
% INPUT
%       V               The metric to represent on y-axis in final plot
%       sec             A structure with fields:
%                       - def
%                       - vec
%                       - col
%                       The field 'vec' may however be omitted but the 
%                       other ones are needed.
%       correct_index   1-by-M cell consisting of correct indices of sites

%%
N_pts=numel(V); N_sec=numel(sec);
% Default values for non-specified fields
if ~isfield(sec,'vec')
    for i=1:N_sec, sec(i).vec=ones(N_pts,1); end
end

map=.8*ones(N_pts,3);

for s=1:N_sec
    zp=sec(s).vec/max(sec(s).vec);
    for p=1:numel(sec(s).def)
        i=sec(s).def(p);
            map(i,:)=sec(s).col;
    end
end

len=numel(correct_index);
hold on;
for i=1:len
    hbar(i)=bar(i,V(i),.8);
    set(hbar(i),'FaceColor',map(i,:),'EdgeColor',map(i,:));
end;
hold off
axis([0 len+1 0 1]);grid on
set(gca,'XTick',[1:30:len]);
set(gca,'XTickLabel',correct_index([1:30:len]));
set(gca,'YTickLabel','');