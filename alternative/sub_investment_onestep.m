function [pol_kp_unc] = sub_investment_onestep(V1,k_grid,pi_x,theta,delta,q,psi)

% We compute the investment policies of unconstrained firms. This requires
% the computation of k_star1 and k_star2.
% OUTPUTS:
%   "pol_kp_unc"   k'(k,x) for unconstrained firms. This is the value (not
%   index, and it is not on the grid k_grid.
%-------------------------------------------------------------------------%

nk = length(k_grid);
nx = size(pi_x,1);

aux1 = -(1-q*psi*theta*(1-delta));
aux2 = -(1-q*psi*(1-delta))*theta;

k_star1 = zeros(nx,1);
k_star2 = zeros(nx,1);
%k_star1_ind = ones(nx,1);
%k_star2_ind = ones(nx,1);
kprime_vec = k_grid;
% TODO: golden maximization?
for x_c = 1:nx
    EVx = zeros(nk,1);
    for xp_c = 1:nx
        EVx = EVx+pi_x(x_c,xp_c)*max(theta*(1-delta)*kprime_vec,V1(:,xp_c));
    end %end x'
    RHS1 = aux1*kprime_vec+q*(1-psi)*EVx;
    RHS2 = aux2*kprime_vec+q*(1-psi)*EVx;
    [~,max_ind1] = max(RHS1);
    [~,max_ind2] = max(RHS2);
    k_star1(x_c) = k_grid(max_ind1);
    k_star2(x_c) = k_grid(max_ind2);
    %k_star1_ind(x_c) = max_ind1;
    %k_star2_ind(x_c) = max_ind2;
end

% Compute k'(k,x), capital investment for unconstrained firms
pol_kp_unc = zeros(nk,nx); % value
for x_c = 1:nx
    for k_c = 1:nk
        k_val = k_grid(k_c);
        if (1-delta)*k_val>k_star2(x_c)
            pol_kp_unc(k_c,x_c) = k_star2(x_c);
        elseif (1-delta)*k_val>=k_star1(x_c) && (1-delta)*k_val<=k_star2(x_c)
            pol_kp_unc(k_c,x_c) = (1-delta)*k_val;
        else
            pol_kp_unc(k_c,x_c) = k_star1(x_c);
        end
    end
end

end %end function "sub_investment_onestep"
