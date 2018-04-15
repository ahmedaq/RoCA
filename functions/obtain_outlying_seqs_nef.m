function outlying_seqs = obtain_outlying_seqs_nef(msa)

% Code to identify the clusters of outlying sequences in HIV Nef
% 
% Last updated: 2017-03-25
% 
% INPUT
%     msa                   Downloaded amino acid MSA
%
% OUTPUT
%     outlying_seqs         Indices of outlying sequences

%% Calculating the similarity matrix

Code_aminoacid='ACDEFGHIKLMNPQRSTVWY';

msa_n = zeros(size(msa));
for kk = 1:length(Code_aminoacid)
    msa_n((msa==Code_aminoacid(kk)))=kk;
end

[Nstar,ls]=size(msa_n); 
msa_extended=zeros(Nstar,20*ls);
for kk=1:ls 
    for Num_aminoacid=1:20 
        msa_extended(:,20*(kk-1)+Num_aminoacid)=(msa_n(:,kk)==Num_aminoacid); 
    end
end
Gamma = (msa_extended*msa_extended')/ls; 

%% Eigenvalue decomposition of the similarity matrix

eig_vec = eig_sort(Gamma);

%% Identifying outlying sequences by visual inspection
outlying_seq1 = find(eig_vec(:,1)>-0.014); 
outlying_seq2 = find(eig_vec(:,2)<-0.04); 
outlying_seqs = [outlying_seq1; outlying_seq2];

% figure;
% scatter(eig_vec(:,1),eig_vec(:,2));
% xlabel('Eigenvector 1')
% ylabel('Eigenvector 2')
% axis([-0.1 0.1 -0.1 0.1])
% grid on
% hold on
% scatter(eig_vec1(outlying_seqs),eig_vec2(outlying_seqs),'r*')
