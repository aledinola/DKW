function [kp_bar] = sub_kp_onestep(profit_mat,b_grid,k_grid,q,theta,delta,lambda)

% Compute the upper bound for k' for constrained firms

[nk,nx] = size(profit_mat);
nb      = size(b_grid,2);

kp_bar   = zeros(nk,nb,nx);

for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            k_val      = k_grid(k_c);
            b_val      = b_grid(k_c,b_c);
            profit_val = profit_mat(k_c,x_c);
            % If k'>=(1-delta)*k
            kp_bar_up   = (profit_val+(1-delta)*k_val-b_val)/(1-q*lambda);
            % If k'<(1-delta)*k
            kp_bar_down = (profit_val+theta*(1-delta)*k_val-b_val)/(theta-q*lambda);
            if profit_val+q*lambda*(1-delta)*k_val-b_val>=0
                kp_bar(k_c,b_c,x_c) = kp_bar_up;
            else
                kp_bar(k_c,b_c,x_c) = kp_bar_down;
            end
        end
    end
end

end %end function "sub_kp_onestep"