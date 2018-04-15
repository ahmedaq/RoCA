function [msa_binary,freq_binary,prev_aa] = binary_approx(msa)

%Code to convert amino acid MSA into binarized MSA
% 
% (adapted from Halabi et al. [1])
% [1] Halabi, N. et al. Supplemental Data Theory Protein Sectors: 
%     Evolutionary Units of Three-Dimensional Structure. Cell 138, (2009).
% 
% INPUT
%       msa             N-by-M amino acid MSA matrix
%
% OUTPUT
%       msa_bin         Binarized MSA with most prevalent amino acid at 
%                       each site represented by 0 and mutant by 1
%       freq_bin        Frequency of most prevalent amino acid at each 
%                       site M

%%
[N,M] = size(msa);

Code_aminoacid='ACDEFGHIKLMNPQRSTVWYB';

%% Frequencies of amino acids at each position
freq=zeros(21,M);
for a=1:21 
    freq(a,:)=sum(msa==Code_aminoacid(a))./N; 
end
% freq(a,i) gives the frequency of amino acid a at position i.

[freq_binary,prev_aa] = max(freq(1:20,:));
% prev_aa is the consensus amino acid at each position and
% freq_bin is its associated frequency

%% Generating binary MSA

msa_binary=1.*(msa~=repmat(Code_aminoacid(prev_aa),N,1));
% msa_bin is the binarized MSA with 0 representing consensus amino acid and
% 1 representing all other possibilities at this position

%% Entropy using binary approximation
Entropy_binary= -freq_binary.*log(freq_binary)-(1-freq_binary).*log(1-freq_binary);

%% Overall entropy
Entropy_complete = zeros(1,M);
for i=1:M
    for a=1:21
        if(freq(a,i)>0)
            Entropy_complete(i)=Entropy_complete(i)-freq(a,i)*log(freq(a,i));
        end
    end
end

%% Generating figure

grey_color = 200*ones(1,3)/255;
figure('name','Validation of binary approximation of MSA'); 
plot(Entropy_complete,Entropy_binary,'o','MarkerFaceColor',grey_color,'MarkerEdgeColor',115*ones(1,3)/255,'MarkerSize',8);
hold on
plot([0:0.001:2],[0:0.001:2],'k')
% figure(3); plot(D_glo,D_bin_wg,'o');
xlabel('$${\rm H}_{i}$$','interpreter','latex');
ylabel('$${\rm H}_{i}^{\rm bin}$$','interpreter','latex');
axis([0 2 0 2])
[pearson_corr_bin_approx,pvalue_bin_approx] = corr(Entropy_complete(:),Entropy_binary(:),'type','Pearson');

fprintf('Pearson correlation of %.2f (P-value = %.1e) is observed between entropy of MSA and that of binarized MSA. \n', ...
    pearson_corr_bin_approx,pvalue_bin_approx);

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...  
  'XColor'      , [0.1 0.1 0.1], ...
  'YColor'      , [0.1 0.1 0.1], ...  
  'YTick'       , 0:.4:2, ...
  'XTick'       , 0:.4:2, ...
  'LineWidth'   , 1        );
