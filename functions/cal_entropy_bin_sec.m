function Entropy_bin_sec = cal_entropy_bin_sec(sec,freq_bin)

% Code to calculate entropy of sector assuming independent sites
%
% Last updated: 2017-03-25
% 
% INPUT
%     sec                   Sector sites
%     freq_bin              Frequency of most prevalent amino acid 
%                           at each site
%
% OUTPUT
%     Entropy_bin_sec       Entropy of sector assuming independent sites

%%

Entropy_bin_sec = sum(- freq_bin(sec).*log(freq_bin(sec))-...
    (1-freq_bin(sec)).*log(1-freq_bin(sec)));
