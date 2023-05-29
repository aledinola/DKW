function [sol,agg,b_grid,distribS,prices,model_mom,flag_ss,par] = fun_steady_state(par)
% fun_steady_state solves the steady-state model given parameter values par
% (including grids,etc.)
% Overview:
%     Step 0: given par, set grids etc.
%     Step 1: given prices (q,w,R), solve firm's dynamic programming problem
%     Step 2: recover firm's value function, entry and exit rule
%     Step 3: calculate stationary distributions mu (mu0 in the draft) and
%     mu_active (mu in the draft)
%     Step 4: Market clearing to find aggregates
%     Step 5: Compute model moments
% Notes:
%     For an explanation of the numerical algorithm to solve for the
%     steady-state, please read Appendix F of our paper.


%------------------------- input checks ----------------------------------%
if isstruct(par)==0
    error('Input par in fun_solve_model must be a structure!')
end
%-------------------------------------------------------------------------%

flag_ss = 0; % so far so good
verbose = par.verbose;

%% Step 0a: Set x-grids, transition matrix for x and stationary distrib for x
if par.x_process == 1
    %-- AR1
    % Auxiliary parameter
    par.mean_x             = (1-par.rhox)*log(par.x0);
    [par.pi_x, log_x_grid] = markovapprox(par.rhox,par.epsx,par.mean_x,par.cover,par.nx);
    par.x_grid             = exp(log_x_grid);

    % Compute stationary distribution of Markov chain x'|x and call it x_prob
    [eig_vectors,eig_values] = eig(par.pi_x');
    [~,arg]                  = min(abs(diag(eig_values)-1));
    unit_eig_vector          = eig_vectors(:,arg);
    par.x_prob               = unit_eig_vector/sum(unit_eig_vector);

    % Compute x-distribution for new entrants (prod shifted by xi, called x0 in the draft)
    [par.x0_prob] = fun_x_entrants(par.x_grid,par.x_prob,par.epsx,par.rhox,par.mean_x,par.xi);

elseif par.x_process == 2
    %-- Bounded Pareto with persistence
    par.x_grid = linspace(par.x_lb,par.x_ub,par.nx)';
    [par.x_prob,par.pi_x] = paretojo(par.nx,par.x_grid',par.x_shape,par.x_rho);
    par.x_prob = par.x_prob';
    % TODO: how to shift the distribution for entrants??
    par.x0_prob = par.x_prob;
else
    error('fun_steady_state: x_process is invalid!')
end

%% Step 0b: Set grid for capital and prob distribution for initial k

% par.k_grid is created in set_parameters.m

% Compute k-distribution for new entrants
if par.k_distrib == 1 % uniform distribution
    par.prob_k = fun_k_entrants_uniform(par.k_grid,par.k_min,par.k_max);
elseif par.k_distrib == 2 % pareto distribution
    par.prob_k = fun_k_entrants_pareto(par.k_grid,par.k_min,par.k_alpha);
else
    error('fun_steady_state: k_distrib is invalid!')
end

%% Step 0c: set auxiliary parameters

par.lambda = par.lambda0*par.theta*(1-par.delta_k);

%% Step 0d: Set vector for fixed cost (depends on kappa)

par.fixcost = zeros(par.nx,1);
for k_c = 1:par.nk
    k_val = par.k_grid(k_c);
    par.fixcost(k_c) = fun.fun_fixcost(k_val,par.fixcost1,par.fixcost2,par);
end

if par.verbose>=1
    fprintf('  \n')
    disp('--------------------------------------------')
    disp('Grid dimensions')
    disp('--------------------------------------------')
    fprintf('nb:   %d  \n',par.nb)
    fprintf('nx:   %d  \n',par.nx)
    fprintf('nk:   %d  \n',par.nk)
    fprintf(' \n')
end

%% Step 1: prices
[prices] = fun_prices(par);

%% Step 2: Value function iter

if verbose>=1; tic; end
[sol,b_grid,phi_dist,flag_vf] = fun_vfi1(prices,par); 
if verbose>=1
    time=toc;
    fprintf('Time to do VFI: %8.4f \n', time)
    fprintf(' \n')
end
if flag_vf<0 
    warning('Some error occurred in fun_vfi1');
    flag_ss = -1;
    sol=[];agg=[];b_grid=[];distribS=[];prices=[];model_mom=[];
    return
end
if any(isnan(phi_dist(:))) || any(isinf(phi_dist(:)))
    warning('phi_dist has NaN/Inf values');
    flag_ss = -1;
    sol=[];agg=[];b_grid=[];distribS=[];prices=[];model_mom=[];
    return
end

%% Distribution
if verbose>=1; tic; end
[mu,mu_active,entry_vec,flag_mu,dist,iter_mu] = fun_distrib1(par,sol,b_grid,phi_dist);

if verbose>=1
    time=toc;
    fprintf('Time to do DISTRIBUTION: %8.4f \n', time)
end
if flag_mu<0
    warning('Distrib did not converge!');

    fprintf('Iter so far = %d    \n',iter_mu)
    fprintf('Last err    = %.15f \n',dist)

    flag_ss = -1;
    %sol=[];agg=[];b_grid=[];distribS=[];prices=[];model_mom=[];
    %return
end
% Group outputs of distribution into a struct:
distribS.mu        = mu;
distribS.mu_active = mu_active;
distribS.entry_vec = entry_vec;

%% Aggregate variables and model moments

if verbose>=1; tic; end
[agg] = fun_aggregates(par,sol,distribS,phi_dist,prices);
if verbose>=1
    time=toc;
    fprintf('Time to do aggregates: %8.4f \n', time)
    fprintf(' \n')
end

% Add phi_dist to <sol> struct
sol.phi_dist = phi_dist;

% Compute model moments (in addition to <agg>)
if verbose>=1; tic; end
model_mom = fun_targets(sol,distribS,par,prices,agg,b_grid);
if verbose>=1
    time=toc;
    fprintf('Time to do targets: %8.4f \n', time)
    fprintf(' \n')
end

if (isstruct(sol)==0)
    error('Output <sol> must be a structure')
end
if (isstruct(agg)==0)
    error('Output <agg> must be a structure')
end
if (isstruct(distribS)==0)
    error('Output <distribS> must be a structure')
end

end %END FUNCTION <fun_steady_state>

