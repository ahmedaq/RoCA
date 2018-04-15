function msa_bin = generate_bin_data_EPmethod(Nbar,C,mu)

% Code to generate Nbar binary data samples satisfying the covariance C 
% and mean mu
% 
% Last updated: 2017-03-25

%%  Generating binary data using the EP method
e = 1;
[msa_bin,~,~,e] = sampleDichGauss01(mu,C,Nbar);   % generate samples from the DG model

if e==0
    fprintf('\n Successfully generated correlated binary samples... \n\n')
else
    fprintf('\nERROR: The specified covariance is NOT FEASIBLE!!')
end

msa_bin = double(msa_bin.');
