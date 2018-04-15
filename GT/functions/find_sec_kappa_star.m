function [sec_kappa_star,kappa_star] = find_sec_kappa_star(num_eig,sec_eig,supp_actual,M)

% Code to find the best alignment (kappa*) of sectors to ground truth units
% 
% Last updated: 2017-03-25

%%

all_possible_perms = perms(1:num_eig);

sec_eig2 = cell(1,size(all_possible_perms,1));
tpr = zeros(size(all_possible_perms,1),num_eig);
fpr = zeros(size(all_possible_perms,1),num_eig);
acc = zeros(size(all_possible_perms,1),num_eig);

for kk = 1:size(all_possible_perms,1)
    perm_iter = all_possible_perms(kk,:);
    for mm = 1:num_eig
        sec_eig2{kk}{mm} = sec_eig{perm_iter(mm)};
        tpr(kk,mm) = sum(ismember(supp_actual{mm},sec_eig2{kk}{mm}))/length(supp_actual{mm});
        fpr(kk,mm) = sum(~ismember(sec_eig2{kk}{mm},supp_actual{mm}))/length(setdiff(1:M,supp_actual{mm}));
        acc(kk,mm) = (sum(ismember(supp_actual{mm},sec_eig2{kk}{mm})) + ...
            sum(ismember(setdiff(1:M,supp_actual{mm}),setdiff(1:M,sec_eig2{kk}{mm}))))/M;
    end
end

[~,kappa_star] = max(mean(acc,2));
sec_kappa_star = sec_eig2{kappa_star};