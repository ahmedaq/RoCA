function sec_ord = order_sectors_according_C(sec,C)

% Code to order sites in the sector according to the sum of individual
% sites correlations with other sites in the sector
% 
% Last updated: 2017-03-25
% 
% INPUT
%     sec                       Sector
%     C                         Correlation matrix
%
% OUTPUT
%     sec_ord                   Ordered sector sites

%%

if isempty(sec)
    sec_ord = [];
else
    sum_corr_k = zeros(1,length(sec));
    for k = 1:length(sec)
%         sum_corr_k(k) = sum((C(sec(k),setdiff(sec,sec(k)))));
        sum_corr_k(k) = sum(abs(C(sec(k),setdiff(sec,sec(k)))));
    end
    [~,b]=sort(sum_corr_k,'descend');
    sec_ord = sec(b);
end