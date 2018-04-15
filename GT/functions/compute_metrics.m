function [mean_tpr,mean_fdr,PM_max] = compute_metrics(sec_kappa_star,supp_actual,num_eig)

% Code to compute the metrics used for comparison of PCA and SPCA methods
% 
% Last updated: 2017-03-25

%% Computing mean TPR and mean FDR

tpr = zeros(1,num_eig);
fdr = zeros(1,num_eig);

for mm = 1:num_eig
    tpr(mm) = sum(ismember(supp_actual{mm},sec_kappa_star{mm}))/length(supp_actual{mm});
    fdr(mm) = sum(~ismember(sec_kappa_star{mm},supp_actual{mm}))/length(sec_kappa_star{mm});
end

mean_tpr = mean(tpr);
fdr(isnan(fdr))=[]; %removing NaN values
mean_fdr = mean(fdr);

%% Computing PM_max

PM = [];
for kk = 1:num_eig
    setdiffkk = setdiff(1:num_eig,kk);
    for mm = 1:length(setdiffkk)
        PM = [PM;sum(ismember(supp_actual{kk},sec_kappa_star{setdiffkk(mm)}))/length(supp_actual{kk})*100];
    end
end

PM_max = max(PM);