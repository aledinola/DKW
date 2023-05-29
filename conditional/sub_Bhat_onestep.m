function [B_hat_new,pol_bp_unc] = sub_Bhat_onestep(B_hat,pol_kp_unc,profit_mat,...
    x_tilde_val,k_grid,x_grid,q,theta,delta,lambda)

% DESCRIPTION:
% This function finds B_hat(k,x), the highest debt level consistent with
% non-negative dividends, by iterating on equations (21) and (22) (fixed
% point iteration). It also returns the borrowing policy b'(k,x).
% INPUTS:
%   "B_hat"        Initial maximum debt,                   dim: (nk,nx)
%   "pol_kp_unc"   k'(k,x) policy for unconstrained firms, dim: (nk,nx)
%   "profit_mat"   Static profit,                          dim: (nk,nx)
%   "x_tilde"      Exit cutoff,                            dim: (nk,1)
%   "k_grid"       Fixed grid for capitaldim,              dim: (nk,1)
%   "q,theta,delta,lamdba"                                 dim: scalars
% OUTPUTS:
%   "B_hat_new"    Updated,                                dim: (nk,nx)
%   "pol_bp_unc"   b'(k,x) policy for unconstr. firms,     dim: (nk,nx)
% NOTES:
% sub_Bhat_onestep is also implemented as MEX in sub_Bhat_onestep_mex
% AUXILIARY:
% It calls the MEX file myinterp1q.
%-------------------------------------------------------------------------%


[nk,nx] = size(B_hat);

k_min = k_grid(1);
k_max = k_grid(nk);

%B_hat_new    = zeros(nk,nx);
pol_bp_unc   = zeros(nk,nx);
%B_hat_interp = zeros(1,nx);

kp_mat = max(min(pol_kp_unc,k_max),k_min);

for x_c=1:nx
    kp_val = kp_mat(:,x_c); %dim: (nk,1)
    B_hat_interp = myinterp1q(k_grid,B_hat,kp_val); %dim: (nk,nx)
    x_cut = myinterp1q(k_grid,x_tilde_val,kp_val); % nk,1
    viable_x = x_grid'>=x_cut; % vector of logicals % (nk,nx)
    B_hat_interp(~viable_x) = nan;
    pol_bp_unc(:,x_c) = min(lambda*kp_val,min(B_hat_interp,[],2));
    
end %x

pol_bp_unc(isnan(pol_bp_unc)) = max(min(lambda*pol_kp_unc(isnan(pol_bp_unc)),k_max),k_min);

B_hat_new = profit_mat+q*pol_bp_unc-fun.adjcost(kp_mat,k_grid,theta,delta);

end %end function "sub_Bhat_onestep"

% % NON VECTORIZED VERSION, coded in the MEX file
% for x_c=1:nx
%     for k_c=1:nk
%         k_val  = k_grid(k_c);
%         kp_val = kp_mat(k_c,x_c); %dim: scalar
%         for xp_c = 1:nx
%             B_hat_interp(xp_c) = myinterp1(k_grid,B_hat(:,xp_c),kp_val);
%         end
%         x_cut = myinterp1(k_grid,x_tilde_val,kp_val); % scalar
%         viable_x = x_grid>=x_cut; % vector of logicals % (nx,1)
%         if (any(viable_x==1))
%           pol_bp_unc(k_c,x_c) = min(lambda*kp_val,min(B_hat_interp(viable_x)));
%         else
%             % There are no x's such that x>=x_cut
%             pol_bp_unc(k_c,x_c) = max(min(lambda*pol_kp_unc(k_c,x_c),k_max),k_min);
%         end
%         B_hat_new(k_c,x_c) = profit_mat(k_c,x_c)+q*pol_bp_unc(k_c,x_c)-fun.adjcost_scal(kp_val,k_val,theta,delta);
%
%     end %k
% end %x
%



