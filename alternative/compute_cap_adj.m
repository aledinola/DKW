function [capadj] = compute_cap_adj(pol_kp,mu,mu_active,phi_dist, ...
    pol_entry,pol_exit,mass,k_grid,delta,psi)
% Description: we compute the four types of capital adjustments
% 1. Upward adjustment by active firms
% 2. Downward adjustment by active firms
% 3. capital bought by entrants
% 4. capital sold by exiters
%   
% Inputs
% - pol_kp: capital next period policy, dim (nk,nb,nx)
% - mu,mu_active,phi_dist: distributions (incumbents, active firms,
% entrants), dim (nk,nb,nx)
% - pol_entry,pol_exit: entry and exit policy, dim (nk,nb,nx)
% - mass: scalar, mass of potential entrants
% - k_grid: dim nk
% - delta: scalar, capital depreciation rate
% - psi: scalar, exogenous exit rate
% 
% Outputs
% - capadj: vector 4 by 1

[~,nb,nx] = size(pol_kp);

k_arr  = repmat(k_grid,1,nb,nx);


up_active = max(pol_kp-(1-delta).*k_arr,0);
down_active = max((1-delta).*k_arr-pol_kp,0);


capadj = zeros(4,1);
% 1. Upward adjustment by active firms
capadj(1) = sum(up_active.*mu_active,'all');

% 2. Downward adjustment by active firms
capadj(2) = sum(down_active.*mu_active,'all');

% 3. capital bought by entrants
capadj(3) = sum(mass.*k_arr.*pol_entry.*phi_dist,'all');

% 4. capital sold by exiters
exit_all = psi+(1-psi).*pol_exit;
capadj(4) = sum(k_arr.*exit_all.*mu,'all');



end % end function plot_cap_adj