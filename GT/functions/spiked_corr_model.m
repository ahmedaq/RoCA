function [C,mu,supp_actual] = spiked_corr_model(M,r,S,ell_0,ell_1,ell_r,max_1pt,min_1pt)

% Code to generate a spiked correlation model using the specified 
% parameters
% 
% Last updated: 2017-03-25

%% Joint support of sectors
joint_supp = sum(S);
parts = [0 S(1) S(1)+S(2) S(1)+S(2)+S(3) S(1)+S(2)+S(3)+S(4) S(1)+S(2)+S(3)+S(4)+S(5)];

%% Generating 1 point frequencies of each site
random_1pts = (max_1pt-min_1pt).*rand(M,1) + min_1pt; 
sigma_site = (1-(1-2*random_1pts).^2)/4;

%% Defining the support location and asscoiated amplitude of each sector
supp_actual = cell(1,r);
for kk = 1:r
    supp_actual{kk} = parts(kk)+1:parts(kk+1);
end

% within the support of each eigenvector, let's specify the ratio of
% positive (and negative) correlations...
per_pos_values = .6;
ratioPos = per_pos_values*ones(1,r); % all eigenvectors the same in this example
ratioNeg = 1-ratioPos;

% Note: these ratios are the number of positive (or neg) correlations over
% the number of all possible correlations within the support of the
% eigenvector ( i.e. ratioPos = num. pos. correlations / s(s-1)/2 ) where s
% is the support size

% Since the number of posibilities is discrete (determined by the number of positive and
% negative entries), ratioPos (and ratioNeg) are in many cases not exactly feasible.
% Thus, we obtain the closest-to-the-specified solution to determine the
% number of positive and negative entries in each eigenvector. This is done
% by the function posneg.m, which returns, in addition to the number of
% pos. and neg. entries in the support, the corresponding actual ratioPos and ratioNeg.  

supp_sizes = parts(2:length(parts)) - parts(1:length(parts)-1);
numPos = zeros(1,r);
numNeg = zeros(1,r);
for kk=1:r
    [numPos(kk),numNeg(kk),ratioPos(kk),ratioNeg(kk)] = posneg(supp_sizes(kk),ratioPos(kk));
end

if joint_supp < r
    disp('ERROR: sparsity is less than rank!!')
    return
end

%% All eigenvalues
ells = ell_1:-(ell_1-ell_r)/(r-1):ell_r;


%% Constructing the eigenvectors with assigned support S
U = zeros(M,r);
for ii=1:r   
    U(supp_actual{ii},ii) = [ones(numPos(ii),1) ; -1*numPos(ii)/numNeg(ii)*ones(numNeg(ii),1)];
    U(:,ii) = U(:,ii)/norm(U(:,ii)); 
end

%% Incorporating phylogeny
U = [1/sqrt(M)*ones(M,1) U];
D = diag([ell_0 ells(1:r)]);

%% Correlation matrix
R = eye(M) - diag(diag(U*D*U')) + U*D*U';

%% Variances (1-point probabilities)
Sigma = diag(sigma_site);

%% 1-point probabilities obtained from the diagonal of C
mu = 1/2*(1-sqrt(1-4*diag(Sigma)));

%% Covariance matrix (prior to normalization)
C = sqrtm((Sigma)) * R * sqrtm((Sigma));