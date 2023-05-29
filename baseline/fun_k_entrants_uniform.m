function [k_prob] = fun_k_entrants_uniform(k_grid,k_min,k_max)
% fun_k_entrants_uniform computes the k-distribution for new entrants
% Inputs: k_grid (grid for k), k_min,k_max.
% k_min is equal to the lower bound for k_grid
% k_max should be below the upper bound of k_grid.
% k_distrib: indicator for the type of distribution. 

%k_min = k_grid(1);
if k_min<=0 || k_max<=0
    error("Inputs <k_min> and <k_max> must be strictly positive!")
end

nk     = length(k_grid);
k_prob = zeros(nk,1);


for k_c = 1:nk
    k_val = k_grid(k_c);
    if k_val>=k_min && k_val <=k_max
        k_prob(k_c) = 1;
    end
end

% k_prob must sum to one
k_prob = k_prob/sum(k_prob);


end %end function "fun_k_entrants_uniform"
