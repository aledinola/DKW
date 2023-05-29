function [exitSS,exit_nogrant,exit_grant_baseline,bin_lids1] = plot_micro_l(sol,distribS,prices,...
    pol_tran_nogrant,weights_nogrant,pol_tran_grant_baseline,weights_grant_baseline,par)
% DESCRIPTION:
%   plot_micro_l generates exit variables by bins of employment l(k,x) 
% INPUTS:
%   "perct"               Vector of percentiles for binning.
%   "sol,distribS,prices" Structures with policies, distribution and prices
%                         in the steady-state
%   "pol_tran_nogrant"    Structures for no grant economy
%   "pol_tran_grant_baseline" Struct for grant baseline economy
%   "par"                 Structure for parameters
% OUTPUTS
%   "exitSS"              Struct with stats by employment bin in the steady state 
%   "exit_nogrant,"       Struct with stats by employment bin in t=1, no grant economy 
%   "exit_grant_baseline" Struct with stats by employment bin in t=1, base grant economy


if ~isstruct(par)
    error('Input par in plot_micro must be a structure')
end

% Outputs:
% exit variables by x for steady state, no grant, and baseline grant
% scenarios.

nb     = par.nb;     % num. gridpoints for debt 
nx     = par.nx;     % num. gridpoints for productivity
nk     = par.nk;     % num. gridpoints for capital small firms 
nn     = par.nn;     % points for dummy variable "impact x grant"
x_grid = par.x_grid; % grid for productivity, "x"
k_grid = par.k_grid; % grid for capital small firms, "kappa"
psi    = par.psi;    % exogenous exit rate

% Steady state objects
pol_exit_ss  = sol.pol_exit;
mu_ss        = distribS.mu;
mu_active_ss = distribS.mu_active;
% transition objects
pol_exit_nogrant        = pol_tran_nogrant.pol_exit;
pol_exit_grant_baseline = pol_tran_grant_baseline.pol_exit;

% Compute l(k,x)
%fun_l(x,wage,k,par)

%mu_mat = squeeze(sum(mu_active_ss,2));
%mu_mat = mu_mat/sum(mu_mat,"all");

l_mat = zeros(nk,nx);
for x_c = 1:nx
    for k_c = 1:nk
        wage = prices.wage;
        k_val = k_grid(k_c);
        x_val = x_grid(x_c);
        l_mat(k_c,x_c) = fun.fun_l(x_val,wage,k_val,par);
    end
end

% Sort by employment
%l_vec  = l_mat(:);
%mu_vec = mu_mat(:);
%perct       = [0.25, 0.5, 0.75]'; 
%bin_lids = quantili(l_vec,mu_vec,perct);
%nbins = length(perct)+1;
%bin_lids1 = [0;bin_lids;inf];

bin_lids1 = [0;10;20;100;inf];
nbins = length(bin_lids1)-1;
bin_mat = max(1,min(discretize(l_mat,bin_lids1),nbins));

%% Steady state
% Use mu0(k,b,x) measure of firms before exit
% Recall: pol_exit has dim: (k,b,x)
exits        = zeros(nbins,1); % measure of exits (without dividing by mu)
emp_share    = zeros(nbins,1); % employment 
output_share = zeros(nbins,1); % output
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            lval = l_mat(k_c,x_c);
            dexit = psi + (1-psi)*pol_exit_ss(k_c,b_c,x_c);
            % Which bin does lval belong to?
            bin_c = bin_mat(k_c,x_c);
            % Update statistics
            % Measure of exits
            exits(bin_c)   = exits(bin_c)+ dexit*mu_ss(k_c,b_c,x_c);
            % Total emp 
            emp_share(bin_c) = emp_share(bin_c) + lval*mu_active_ss(k_c,b_c,x_c);  
            % Total output
            fc = 0;%par.fixcost(k_c);
            output_share(bin_c) = output_share(bin_c) + fun.prod_small(x_grid(x_c),k_grid(k_c),lval,fc,par)* ...
                mu_active_ss(k_c,b_c,x_c);
        end
    end
end
% Share of exits, emp, and output for each bin of emp, dim: (nbins,1)
exitSS.exit_share   = exits/sum(exits); 
exitSS.emp_share    = emp_share/sum(emp_share); 
exitSS.output_share = output_share/sum(output_share); 
exitSS.exits        = exits;

%% Transition, consider only t=1
% measure of exits (without dividing by mu)
exits_nogrant        = zeros(nbins,1); 
mass_nogrant         = zeros(nbins,1);
exits_grant_baseline = zeros(nbins,1); 
mass_grant_baseline  = zeros(nbins,1);

for k_c = 1:nk
    for n_c = 1:nn
        for x_c = 1:nx
            for b_c = 1:nb % current debt
                % Which bin does the firm belong to?
                bin_c = bin_mat(k_c,x_c);
                
                % No grant
                dexit = psi+(1-psi)*pol_exit_nogrant(k_c,b_c,x_c,1,n_c);

                % Measure of exits (recall that the incumbent distribution
                % at t=1 is the same as the steady state incumbent
                % distribution.)
                exits_nogrant(bin_c) = exits_nogrant(bin_c)+dexit*weights_nogrant(k_c,x_c,n_c)*mu_ss(k_c,b_c,x_c);

                % Mass of firm in each bin
                mass_nogrant(bin_c) = mass_nogrant(bin_c) + ...
                    weights_nogrant(k_c,x_c,n_c)*mu_ss(k_c,b_c,x_c);

                % Baseline grant
                dexit = psi+(1-psi)*pol_exit_grant_baseline(k_c,b_c,x_c,1,n_c);

                % Measure of exits (recall that the incumbent distribution
                % at t=1 is the same as the steady state incumbent
                % distribution.)
                exits_grant_baseline(bin_c) = exits_grant_baseline(bin_c)+ ...
                    dexit*weights_grant_baseline(k_c,x_c,n_c)*mu_ss(k_c,b_c,x_c);

                % Mass of firm in each bin
                mass_grant_baseline(bin_c) = mass_grant_baseline(bin_c) + ...
                    weights_grant_baseline(k_c,x_c,n_c)*mu_ss(k_c,b_c,x_c);
            end % b_c
        end % x_c
    end % n_c
end % k_c

exit_nogrant.exit_share = exits_nogrant/sum(exits_nogrant); 
exit_nogrant.exits      = exits_nogrant;
exit_nogrant.exit_rate  = exits_nogrant./mass_nogrant;

exit_grant_baseline.exit_share = exits_grant_baseline/sum(exits_grant_baseline); 
exit_grant_baseline.exits      = exits_grant_baseline;
exit_grant_baseline.exit_rate  = exits_grant_baseline./mass_grant_baseline;

end %END FUNCTION "plot_micro_l"

