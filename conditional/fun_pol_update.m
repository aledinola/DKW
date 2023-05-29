function [pol_debt,pol_kp,pol_kp_ind,val] = fun_pol_update(val,val_unc,pol_bp_unc,...
    pol_kp_unc,pol_kp_ind_con,profit_mat,B_hat,k_grid,b_grid,q,theta,delta)

[nk,nb,nx] = size(val);

% Compute optimal debt policy implied by eq. (25)
pol_debt   = zeros(nk,nb,nx);
pol_kp     = zeros(nk,nb,nx);
pol_kp_ind = ones(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            k_val      = k_grid(k_c);
            b_val      = b_grid(k_c,b_c);
            profit_val = profit_mat(k_c,x_c);
            if b_val<=B_hat(k_c,x_c) % firm is unconstrained
                val(k_c,b_c,x_c) = val_unc(k_c,b_c,x_c);
                pol_debt(k_c,b_c,x_c) = pol_bp_unc(k_c,x_c);
                pol_kp(k_c,b_c,x_c) = pol_kp_unc(k_c,x_c);
                % Left grid point
                %pol_kp_ind(k_c,b_c,x_c) = find_loc(k_grid,pol_kp(k_c,b_c,x_c));
                % Closest grid point
                [~,pol_kp_ind(k_c,b_c,x_c)] = min(abs(k_grid-pol_kp(k_c,b_c,x_c)));
                %pol_kp_ind(k_c,b_c,x_c) = pol_kp_ind_con(k_c,b_c,x_c);
            else % firm is constrained
                
                kprime = k_grid(pol_kp_ind_con(k_c,b_c,x_c));
                bprime = max(b_grid(pol_kp_ind_con(k_c,b_c,x_c),1),...
                    (1/q)*(b_val-profit_val+fun.adjcost_scal(kprime,k_val,theta,delta)));
                %bprime = (1/q)*(b_val-profit_val+fun.adjcost_scal(kprime,k_val,theta,delta));
                pol_debt(k_c,b_c,x_c) = bprime;
                pol_kp(k_c,b_c,x_c)   = kprime;
                pol_kp_ind(k_c,b_c,x_c) = pol_kp_ind_con(k_c,b_c,x_c);
            end
        end
    end
end

end %end function "fun_pol_update"