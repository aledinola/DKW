function [pol_entry,pol_exit,pol_exit_forced,pol_exit_vol] = ...
    fun_entry_exit(val,profit_mat,b_grid,k_grid,theta,delta,cost_e)

% DESCRIPTION
% Compute entry and exit policy functions on (k,b,x)

[nk,nb,nx] = size(val);

% Entry policy d^e(k,b,x) 
pol_entry  = zeros(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            k_val = k_grid(k_c);
            b_val = b_grid(k_c,b_c); 
            aux1 = val(k_c,b_c,x_c)>=cost_e+k_val-b_val;
            %aux2 = profit_vec(x_c,k_c)-b_val+q*theta*kappa>=0;
            % Maybe aux2 is important!!!
            %aux = (aux1 && aux2);
            pol_entry(k_c,b_c,x_c) = double(aux1);
        end
    end
end

% Exit policy d^l(x,b)
pol_exit        = zeros(nk,nb,nx);
pol_exit_vol    = zeros(nk,nb,nx);
pol_exit_forced = zeros(nk,nb,nx);

for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            kappa = k_grid(k_c);
            b_val = b_grid(k_c,b_c);
            % Forced liquidation
            aux1 = profit_mat(k_c,x_c)-b_val+theta*(1-delta)*kappa<0;
            pol_exit_forced(k_c,b_c,x_c) = double(aux1);
            % Voluntary liquidation
            aux2 = val(k_c,b_c,x_c)< theta*(1-delta)*kappa-b_val;
            pol_exit_vol(k_c,b_c,x_c) = double(aux2);
            aux  = (aux1 || aux2);
            pol_exit(k_c,b_c,x_c) = double(aux);
        end
    end
end


end %end function "fun_entry_exit"