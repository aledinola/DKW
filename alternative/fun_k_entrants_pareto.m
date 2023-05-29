function [k_prob] = fun_k_entrants_pareto(k_grid,k_min,k_alpha)
% This function computes the k-distribution for new entrants 
% Inputs: k_grid (grid for k), k_alpha (curvature parameter of the pareto
% distribution.

% k_alpha and k_min must be positive
 
%k_min = k_grid(1);
if k_alpha<=0 || k_min<=0
    error("Input <k_alpha> and <k_min> must be strictly positive!")
end

nk     = length(k_grid);
k_prob = zeros(nk,1);


for k_c = 1:nk
    k_val = k_grid(k_c);
    if k_val>=k_min
        % This is the p.d.f. of the Pareto distribution
        k_prob(k_c) = k_alpha*k_min^k_alpha/k_val^(1+k_alpha);
    end
end

% k_prob must sum to one
k_prob = k_prob/sum(k_prob);


end %end function