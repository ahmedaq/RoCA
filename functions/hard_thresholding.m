function [value_after_thresh] = hard_thresholding(x,t)

% Zongming Ma's [1] code available at 
% http://www-stat.wharton.upenn.edu/~zongming/software/SPCALab/SPCALab.zip
% [1] Ma, Z. (2013) Sparse principal component analysis and iterative 
% thresholding. Ann Statist, Vol.41, pp.772-801. 

n = length(x);
value_after_thresh = zeros(1,n);
for i = 1:n 
	if x(i) >= t
	  value_after_thresh(i) = x(i);
	elseif x(i) <= -t
	  value_after_thresh(i) = x(i);
	else value_after_thresh(i) = 0;
	end
end