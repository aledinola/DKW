function [entry,output_small,liq,L_small,cost_adj,exit_rate,entry_rate] = ...
    sub_aggregates_onestep(mu,mu_active,pol_kp_ind,pol_entry,exit_all,...
    phi_dist,wage,mass,A,par)

% This subfunction is called by fun_aggregates and fun_aggregates_tran.
% NOTE: A is a scalar and is equal to the productivity shock. In the
% steady-state A=1 by construction.

% Unpack 
k_grid  = par.k_grid;
x_grid  = par.x_grid;
fixcost = par.fixcost;
cost_e  = par.cost_e;
theta   = par.theta;
delta_k = par.delta_k;
nb      = par.nb;
nx      = par.nx;
nk      = par.nk;

%entry        = 0;
%output_small = 0; %output (small firms)
%liq          = 0; %liquidation (small firms)
%L_small      = 0; %employment (small firms)
%cost_adj     = 0; %capital adjustment cost (small firms)
% for x_c = 1:nx
%     x_val = x_grid(x_c);
%     for b_c = 1:nb
%         for k_c = 1:nk
%             kappa        = k_grid(k_c);
%             kappa_next   = k_grid(pol_kp_ind(k_c,b_c,x_c));
%             c            = fixcost(k_c);      
%             l_opt        = fun.fun_l(x_val,wage,kappa,par);
%             % Update L_small, output_small, cost_adj,liq,entry
%             L_small      = L_small + l_opt*mu_active(k_c,b_c,x_c);
%             output_small = output_small + fun.prod_small(x_val,kappa,l_opt,c,par)*mu_active(k_c,b_c,x_c);
%             cost_adj     = cost_adj + fun.adjcost_scal(kappa_next,kappa,theta,delta_k)*mu_active(k_c,b_c,x_c);
%             liq          = liq + theta*(1-delta_k)*kappa*exit_all(k_c,b_c,x_c)*mu(k_c,b_c,x_c);
%             entry        = entry  + mass*(kappa+cost_e)*pol_entry(k_c,b_c,x_c)*phi_dist(k_c,b_c,x_c);
%         end
%     end
% end


kappa = repmat(k_grid,1,nb,nx);
c     = repmat(fixcost,1,nb,nx);
x_val = permute(repmat(A*x_grid,1,nk,nb),[2,3,1]);

l_opt   = fun.fun_l(x_val,wage,kappa,par);         % (nk,nb,nx)
y_opt   = fun.prod_small(x_val,kappa,l_opt,c,par); % (nk,nb,nx)

L_small      = sum(l_opt.*mu_active,'all'); % scalar
output_small = sum(y_opt.*mu_active,'all'); % scalar
liq          = sum(theta*(1-delta_k)*kappa.*exit_all.*mu,'all'); % scalar
% entry is the total cost of entry
entry        = sum(mass*(kappa+cost_e).*pol_entry.*phi_dist,'all'); % scalar
% entry_rate is the measure of entrants
entry_rate   = sum(mass*pol_entry.*phi_dist,'all'); % scalar
% exit_rate is the measure of firms that exit
exit_rate    = sum(exit_all.*mu,'all');  % scalar

cost_adj     = 0; %capital adjustment cost (small firms)

for x_c = 1:nx
    for b_c = 1:nb
        k_val    = k_grid; % (nk,1)
        kp_val   = k_grid(pol_kp_ind(:,b_c,x_c));
        temp     = fun.adjcost(kp_val,k_val,theta,delta_k).*mu_active(:,b_c,x_c);
        cost_adj = cost_adj + sum(temp);
    end
end


end %end function "sub_aggregates_onestep"







