function [sim] = fun_simulate_agesize(sol,prices,par,distrib,b_grid,N_sim,T_sim)
% fun_simulate_agesize simulates a panel of firms from steady-state of the model.
% This function is called from TODO: check main_tables.m
% TODO: Now we simulate only unconstrained firms (we draw the initial conditions
% for (k,b,x) from the stationary distribution but considering only the
% (k,b,x) that are unconstrained. But we want the size-age relationship for
% all firms, also the constrained ones.
% INPUTS:
% sol,prices,par,distrib: steady state structures
% b_grid is 2-dim array (nk,nb)
% N_sim, T_sim are scalars
% OUTPUTS:
% sim is a structure with fields:
% - k_sim, b_sim, x_sim: simulated indexes
% - k_sim_val,x_sim_val: simulated values
% b_sim_val is missing
% AUXILIARY FUNCTIONS:
% - markov_sim
% - find_loc

%% Unpack
nk     = par.nk;
nb     = par.nb;
%nx    = par.nx;
pi_x   = par.pi_x; % Transition matrix exogenous shock (x,x')
k_grid = par.k_grid; % (nk,1)
x_grid = par.x_grid; % (nx,1)

pol_kp_ind = sol.pol_kp_ind; %integers, (k,b,x)
%pol_kp     = sol.pol_kp;
pol_debt   = sol.pol_debt; % (k,b,x)
pol_exit   = sol.pol_exit; % (k,b,x)
pol_entry  = sol.pol_entry; % (k,b,x)
wage       = prices.wage; 
mu_active  = distrib.mu_active; % (k,b,x)
phi_dist   = sol.phi_dist; %(k,b,x)

%% Reset seed
rng(par.seed);

%% Simulate exogenous stochastic productivity x from transition matrix


% Initial distribution for x: from mu_active, but only unconstrained firms
prob0V    = squeeze(sum(phi_dist.*pol_entry,[1,2]));
prob0V    = prob0V/sum(prob0V);
unif_rand = rand(N_sim, T_sim);

% x_sim is (N,T) a panel of integers {1,2,..,NX}
x_sim     = markov_sim(N_sim, T_sim, prob0V, pi_x', unif_rand, 1);
x_sim_val = x_grid(x_sim);

%% Simulate endogenous variables k, b and survival
k_sim     = ones(N_sim,T_sim); % indexes
b_sim     = ones(N_sim,T_sim); % indexes
% Indicator for survival till the end. It is initialized at 1 and switches
% to 0 as soon as exit is true
surv_sim  = ones(N_sim,T_sim);
l_sim  = zeros(N_sim,T_sim); % size

% Draw initial k and b. 
kb_grid = (1:nk*nb)';
for i_c =1:N_sim
    init_mu = phi_dist(:,:,x_sim(i_c,1)).*pol_entry(:,:,x_sim(i_c,1));
    kb_prob = init_mu(:)/sum(init_mu(:));
    % Output of simulate_iid is an index
    kb_ind = simulate_iid_fast(kb_grid,kb_prob,0);
    [k_init,b_init] = ind2sub([nk,nb],kb_ind);
    k_sim(i_c,1) = k_init;
    b_sim(i_c,1) = b_init;
    % Compute firm size in t_c=1
    l_sim(i_c,1) = fun.fun_l(x_grid(x_sim(i_c,1)),wage,k_grid(k_sim(i_c,1)),par);
end

% Simulate for t>1
unif_rand1 = rand(N_sim, T_sim);
for i_c = 1:N_sim
    for t_c = 1:T_sim-1

        if surv_sim(i_c,t_c) ==1 % if firm i_c survives in t_c
            k_sim(i_c,t_c+1) = pol_kp_ind(k_sim(i_c,t_c),b_sim(i_c,t_c),x_sim(i_c,t_c));
            % b_next is a value, does not lie necessarily on the grid
            b_next           = pol_debt(k_sim(i_c,t_c),b_sim(i_c,t_c),x_sim(i_c,t_c));
            b_grid_temp      = b_grid(k_sim(i_c,t_c+1),:)';
            % omega is the weight on the left grid point
            [left_loc,omega] = find_loc(b_grid_temp,b_next);
            u = unif_rand1(i_c,t_c);
            % This is the index for b_next
            b_sim(i_c,t_c+1) = left_loc+(u>omega);
            % Compute firm size
            l_sim(i_c,t_c+1) = fun.fun_l(x_grid(x_sim(i_c,t_c+1)),wage,k_grid(k_sim(i_c,t_c+1)),par);

            %b_sim(i_c,t_c+1) = closest;
            surv_sim(i_c,t_c+1:end) = surv_sim(i_c,t_c)*(pol_exit(k_sim(i_c,t_c),b_sim(i_c,t_c),x_sim(i_c,t_c))==0);
        else
            % Once a firm exits, no need to go on until t=T_sim
            break
        end
    end % for t_c
end % for i_c

k_sim_val = k_grid(k_sim);

% Pack outputs
sim = v2struct(k_sim,b_sim,x_sim,surv_sim,k_sim_val,x_sim_val,l_sim);

end %end function "fun_simulate_agesize"
