function [pol_entry,pol_exit] = interp_entry_exit(pol_entry,pol_exit, ...
    b_grid,val,profit_vec,par)

% See Leo's note "Interpolation of entry and exit policies"

nk    = par.nk;
nb    = par.nb;
nx    = par.nx;
theta = par.theta;
cost_e = par.cost_e;
k_grid = par.k_grid;
delta  = par.delta_k;
entind     = zeros(nk,nb); % for every k,b, finds highest productivity x without entry
valentlow  = zeros(nk,nb); % v_L, net entry value at i (no entry, negative)
valenthigh = zeros(nk,nb); % v_H, net entry value at i+1 (entry, weakly positive)

% entry index
for b_c = 1:nb
    for k_c = 1:nk
        kappa = k_grid(k_c);
        x_i = find(pol_entry(k_c,b_c,:)==0, 1, 'last' ); % first x_index where pol_entry becomes 1 (entry)
        if ~isempty(x_i) && x_i<nx
            entind(k_c,b_c) = x_i;
            valentlow(k_c,b_c)  = val(k_c,b_c,x_i)-cost_e-kappa+b_grid(k_c,b_c);
            valenthigh(k_c,b_c) = val(k_c,b_c,x_i+1)-cost_e-kappa+b_grid(k_c,b_c);
        else
            entind(k_c,b_c) = nx;
        end
    end
end

% exit index
exind = zeros(nk,nb); % for every k and b, finds highest productivity x without entry
valexlow = zeros(nk,nb); % v_L, net entry value at i (no entry, negative)
valexhigh = zeros(nk,nb); % v_H, net entry value at i+1 (entry, weakly positive)

for b_c = 1:nb
    for k_c = 1:nk
        kappa = k_grid(k_c);
        % Find marginal exit index i and net values vL and vH at i and i+1
        x_i = find(pol_exit(k_c,b_c,:)==1, 1, 'last' );
        if ~isempty(x_i) && x_i<nx
            exind(k_c,b_c) = x_i;


            % Forced liq is pivotal
            if profit_vec(k_c,x_i)-b_grid(k_c,b_c)+theta*(1-delta)*kappa<0 ...
                    && profit_vec(k_c,x_i+1)-b_grid(k_c,b_c)+theta*(1-delta)*kappa>=0
                valexlow(k_c,b_c)  = profit_vec(k_c,x_i)-b_grid(k_c,b_c)+theta*(1-delta)*kappa;
                valexhigh(k_c,b_c) = profit_vec(k_c,x_i+1)-b_grid(k_c,b_c)+theta*(1-delta)*kappa;
            elseif val(k_c,b_c,x_i)-theta*(1-delta)*kappa+b_grid(k_c,b_c)<0 ...
                    && val(k_c,b_c,x_i+1)-theta*(1-delta)*kappa+b_grid(k_c,b_c)>=0
                % voluntary liquidation is pivotal
                valexlow(k_c,b_c)  = val(k_c,b_c,x_i)-theta*(1-delta)*kappa+b_grid(k_c,b_c);
                valexhigh(k_c,b_c) = val(k_c,b_c,x_i+1)-theta*(1-delta)*kappa+b_grid(k_c,b_c);
            else
                disp('Neither exit margin relevant')
                keyboard
            end
        else
            exind(k_c,b_c) = nx;
        end
    end
end

% Adjust entry/exit probabilities at margin x_i OR x_i+1
for b_c = 1:nb
    for k_c = 1:nk
        if entind(k_c,b_c)<nx
            sumv = valentlow(k_c,b_c)+valenthigh(k_c,b_c);

            if sumv>0
                pol_entry(k_c,b_c,entind(k_c,b_c)) = sumv/(valenthigh(k_c,b_c)-valentlow(k_c,b_c))/2;
            else
                pol_entry(k_c,b_c,entind(k_c,b_c)+1) = 1-sumv/(valentlow(k_c,b_c)-valenthigh(k_c,b_c))/2;
            end
        end
        if exind(k_c,b_c)<nx
            sumv = valexlow(k_c,b_c)+valexhigh(k_c,b_c);
            if sumv>0
                pol_exit(k_c,b_c,exind(k_c,b_c)) = 1-sumv/(valexhigh(k_c,b_c)-valexlow(k_c,b_c))/2;
            else
                pol_exit(k_c,b_c,exind(k_c,b_c)+1) = sumv/(valexlow(k_c,b_c)-valexhigh(k_c,b_c))/2;
            end
        end
    end
end

end % END function <interp_entry_exit>

