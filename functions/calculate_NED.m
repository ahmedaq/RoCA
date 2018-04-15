function [NED_abs_random,NED_abs_inter,NED_abs_intra] = calculate_NED(msa_bin,PC,alpha)

% Code to calculate Normalized Entropy Deviation (NED) for characterizing
% statistical independence between pairs of inferred sectors
% 
% Last updated: 2017-03-25
% 
% INPUT
%     msa_bin               N-by-M binarized MSA
%     PC                    M-by-(alpha) matrix, columns corresponding to sparse 
%                           estimate of leading eigenvectors
%     alpha                 number of significant eigenvectors
%     
% OUTPUT
%     NED_abs_random        NED for the randomized (null) case 
%     NED_abs_inter         Inter-sector NED 
%     NED_abs_intra         Intra-sector NED

%% Parameters

[N,M] = size(msa_bin);
x_top_elements = 5;
ITER = 500;
no_of_sec_pairs = nchoosek(alpha,2);


%% Inter-sector entropies

fprintf('\n-----------------------------------------------------------------------------------\n')
fprintf('NED calculations\n')
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('Calculating inter-sector entropies...\n')
Entropy_bin_group_sec1 = zeros(alpha-1,alpha);
Entropy_bin_group_sec2 = zeros(alpha-1,alpha);
Entropy_bin_group_sec1_sec2 = zeros(alpha-1,alpha);
Sum_entropy_bin_groups_sec1_sec2 = zeros(alpha-1,alpha);

for kk = 1:alpha-1
    
    [~,b1] = sort(round2dp(abs(PC(:,kk)),5),'descend');
    sec1 = b1(1:x_top_elements);
    
    for mm = kk+1:alpha
        
        [~,b2] = sort(round2dp(abs(PC(:,mm)),5),'descend');
        sec2 = b2(1:x_top_elements);
        
        Entropy_bin_group_sec1(kk,mm) = calc_entropy_bin_group_sec(sec1,msa_bin);
        Entropy_bin_group_sec2(kk,mm) = calc_entropy_bin_group_sec(sec2,msa_bin);
        
        Entropy_bin_group_sec1_sec2(kk,mm) = calc_entropy_bin_group_sec(unique([sec1; sec2]),msa_bin);
        Sum_entropy_bin_groups_sec1_sec2(kk,mm) = Entropy_bin_group_sec1(kk,mm) + Entropy_bin_group_sec2(kk,mm);
        
    end
end

vec_Entropy_group_sec1 = obtain_vec_triu_matrix(Entropy_bin_group_sec1);
vec_Entropy_group_sec2 = obtain_vec_triu_matrix(Entropy_bin_group_sec2);
vec_Entropy_group_sec1_sec2 = obtain_vec_triu_matrix(Entropy_bin_group_sec1_sec2);
sum_vec_Entropy_group_sec1_sec2 = vec_Entropy_group_sec1 + vec_Entropy_group_sec2;


NED_abs_inter = zeros(1,no_of_sec_pairs);
for kk = 1:no_of_sec_pairs
    
    NED_abs_inter(:,kk) = (sum_vec_Entropy_group_sec1_sec2(:,kk)-vec_Entropy_group_sec1_sec2(:,kk))./...
        (vec_Entropy_group_sec1_sec2(:,kk));
end

%% Intra sector entropies

fprintf('Calculating intra-sector entropies...\n')
Entropy_bin_group_sec1_intra = zeros(1,alpha);
Entropy_bin_group_sec2_intra = zeros(1,alpha);
Entropy_bin_group_sec1_sec2_intra = zeros(1,alpha);
Sum_entropy_bin_groups_sec1_sec2_intra = zeros(1,alpha);

for kk = 1:alpha
    
    [~,b1] = sort(abs(PC(:,kk)),'descend');
    sec11 = b1(1:2*x_top_elements);
    sec1 = sec11((1:x_top_elements));
    sec2 = sec11((x_top_elements+1:2*x_top_elements));
    
    Entropy_bin_group_sec1_intra(kk) = calc_entropy_bin_group_sec(sec1,msa_bin);
    Entropy_bin_group_sec2_intra(kk) = calc_entropy_bin_group_sec(sec2,msa_bin);
    
    Entropy_bin_group_sec1_sec2_intra(kk) = calc_entropy_bin_group_sec(([sec1;sec2].'),msa_bin);
    Sum_entropy_bin_groups_sec1_sec2_intra(kk) = Entropy_bin_group_sec1_intra(kk) + Entropy_bin_group_sec2_intra(kk);
    
end

NED_abs_intra = (Sum_entropy_bin_groups_sec1_sec2_intra-Entropy_bin_group_sec1_sec2_intra)./Entropy_bin_group_sec1_sec2_intra;
% nd_sq_intra = (Sum_entropy_bin_groups_sec1_sec2_intra-Entropy_bin_group_sec1_sec2_intra).^2./(Entropy_bin_group_sec1_sec2_intra).^2

%% Computing the NED-Random


fprintf('Calculating randomized case entropies...\n')
NED_abs_random = [];
sum_Entropy_bin_group_sec1_sec2_random = zeros(1,ITER);
Entropy_bin_group_sec1_sec2_random = zeros(1,ITER);

for kk = 1:alpha-1
    for mm = kk+1:alpha
        [~,b1] = sort(abs(PC(:,kk)),'descend');
        sec1_random = b1(1:x_top_elements);
        
        [~,b2] = sort(abs(PC(:,mm)),'descend');
        sec2_random = b2(1:x_top_elements);
        
        parfor iter = 1:ITER
            
            msa_bin_rnd=zeros(N,M);
            perm_seq1=randperm(N);
            msa_bin_rnd(:,sec1_random)=msa_bin(perm_seq1(:),sec1_random);
            perm_seq2=randperm(N);
            msa_bin_rnd(:,sec2_random)=msa_bin(perm_seq2(:),sec2_random);
            
            Entropy_bin_group_sec1_random = calc_entropy_bin_group_sec(sec1_random,msa_bin_rnd);
            Entropy_bin_group_sec2_random = calc_entropy_bin_group_sec(sec2_random,msa_bin_rnd);
            sum_Entropy_bin_group_sec1_sec2_random(iter) = Entropy_bin_group_sec1_random+Entropy_bin_group_sec2_random;
            Entropy_bin_group_sec1_sec2_random(iter) = calc_entropy_bin_group_sec(unique([sec1_random; sec2_random]),msa_bin_rnd);
            
        end
        
        NED_abs_random = [NED_abs_random ...
            ((sum_Entropy_bin_group_sec1_sec2_random-Entropy_bin_group_sec1_sec2_random)./...
            (Entropy_bin_group_sec1_sec2_random)).'];
    end
    fprintf('Pairs of sector %d done...\n',kk)
end

fprintf('Calculation completed\n')