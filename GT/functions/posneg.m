function [P,N,ratioP_out,ratioN_out] = posneg(p,ratioP)

% Code to generate support of eigenvector such that there are 'ratioP'
% positively correlated sites and rest are negatively correlated with these
% 

%%
ratioN = 1-ratioP;

if ((ratioP*(p-1)+1)^2-ratioN^2*(p-1)^2) < 0
    fprintf('\n ERROR: The specified ratio of positive correlations is not feasible!!!\n')
    disp(' Setting ratioP = 1 instead...')
    P = p;
    N = 0;
    ratioP_out = 1;
    ratioN_out = 0;
    return;
end

N = sqrt(p/2*(ratioP*(p-1)+1-sqrt((ratioP*(p-1)+1)^2-ratioN^2*(p-1)^2)));
P = p-N;

P = round(P);
N = p-P;

ratioP_out = (P*(P-1)+N*(N-1))/p/(p-1);
ratioN_out = 2*P*N/p/(p-1);

end

