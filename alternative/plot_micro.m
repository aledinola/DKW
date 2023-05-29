function [exitSS,exit_nogrant,exit_grant_baseline] = plot_micro(sol,distribS,prices,...
    pol_tran_nogrant,distrib_tran_nogrant,path_nogrant,weights_nogrant,...
    pol_tran_grant_baseline,distrib_tran_grant_baseline,path_grant_baseline,weights_grant_baseline,...
    par)

% INPUTS:
%   "sol,distribS,prices" Structures with policies, distribution and prices
%   in the steady-state
%   "pol_tran_nogrant,distrib_tran_nogrant,path_nogrant" Structures for no
%   grant economy
%   "pol_tran_grant_baseline,distrib_tran_grant_baseline,path_grant_baseline"
%   Struct for grant baseline economy
%   "par"

if ~isstruct(par)
    error('Input par in plot_micro must be a structure')
end

% Outputs:
% exit variables by x for steady state, no grant, and baseline grant
% scenarios.

nb     = par.nb;     % num. gridpoints for debt 
nx     = par.nx;     % num. gridpoints for productivity
nk     = par.nk;     % num. gridpoints for capital small firms 
ni     = par.ni;     % points for dummy variable "impact"
ns     = par.ns;     % points for dummy variable "grant"
nn     = par.nn;     % points for dummy variable "impact x grant"
x_grid = par.x_grid; % grid for productivity, "x"
k_grid = par.k_grid; % grid for capital small firms, "kappa"
psi    = par.psi;    % exogenous exit rate
%% Steady-state

% Exit rate as a function of x in the steady state
% Use mu0(k,b,x) measure of firms before exit
% Recall: pol_exit has dim: (k,b,x)
exits_x        = zeros(nx,1); % measure of exits by x (without dividing by mu)
exit_rate_x    = zeros(nx,1); % exits/mass of firms equals exit rate
emp_loss_x     = zeros(nx,1); % employment loss due to exit
emp_share_x    = zeros(nx,1); % employment by x
output_share_x = zeros(nx,1); % output by x
mass_x         = zeros(nx,1); % mass of incumbent firms by x

for x_c = 1:nx
    % Initialize running sums
    exits_x(x_c) = 0;
    emp_loss_x(x_c)=0;
    emp_share_x(x_c)=0;
    output_share_x(x_c)=0;
    for b_c = 1:nb % current debt
        for k_c = 1:nk
            dexit = psi + (1-psi)*sol.pol_exit(k_c,b_c,x_c);
            exits_x(x_c) = exits_x(x_c)+dexit*distribS.mu(k_c,b_c,x_c);
            %             fun_l(x,wage,kappa,par)
            l_small = fun.fun_l(x_grid(x_c),prices.wage,k_grid(k_c),par);
            emp_loss_x(x_c) = emp_loss_x(x_c) + l_small * ...
                dexit*distribS.mu(k_c,b_c,x_c);
            emp_share_x(x_c) = emp_share_x(x_c) + l_small*...
                distribS.mu_active(k_c,b_c,x_c);   %prod_small(x,kappa,l_small,c,par)
            % We set the fixed cost to zero
            fc = 0;%par.fixcost(k_c);
            output_share_x(x_c) = output_share_x(x_c) + fun.prod_small(x_grid(x_c),k_grid(k_c),l_small,fc,par)* ...
                distribS.mu_active(k_c,b_c,x_c);
            mass_x(x_c) = mass_x(x_c) + distribS.mu_active(k_c,b_c,x_c);
        end %k_c
    end %b_c
    exit_rate_x(x_c) = exits_x(x_c)/sum(distribS.mu(:,:,x_c),'all');
end %x_c

% Share of exit due to each x 
exit_share_x   = exits_x/sum(exits_x); %dim: (nx,1)
emp_share_x    = emp_share_x/sum(emp_share_x); %dim: (nx,1)
output_share_x = output_share_x/sum(output_share_x); %dim: (nx,1)

% Save to structure exitSS
exitSS.exit_rate_x  = exit_rate_x;
exitSS.exits_x      = exits_x;
exitSS.exit_share_x = exit_share_x;
exitSS.emp_loss_x   = emp_loss_x;
exitSS.emp_share_x  = emp_share_x;
exitSS.output_share_x = output_share_x;
exitSS.mass_x       = mass_x;


%% Transition without grant, consider only t=1

% State-space: (k,b,x,time,impact x grant)

exits_x     = zeros(nx,1); % measure of exits by x (without dividing by mu)
exit_rate_x = zeros(nx,1); % exits/mass of firms equals exit rate
emp_loss_x  = zeros(nx,1); % employment loss due to exit
mass_x      = zeros(nx,1); % mass of incumbent firms by x
for x_c = 1:nx
    for n_c = 1:nn
        [i_c,~] = ind2sub([ni,ns],n_c); % i_c = impact indicator; s_c = grant indicator
        for k_c = 1:nk
            for b_c = 1:nb % current debt
                % dim: (b,x,kappa,time,impact x grant)
                dexit = psi+(1-psi)*pol_tran_nogrant.pol_exit(k_c,b_c,x_c,1,n_c);
                mass_x(x_c) = mass_x(x_c) + weights_nogrant(k_c,x_c,n_c)*distrib_tran_nogrant.mu(k_c,b_c,x_c,1,n_c);

                exits_x(x_c) = exits_x(x_c)+dexit*weights_nogrant(k_c,x_c,n_c)*distrib_tran_nogrant.mu(k_c,b_c,x_c,1,n_c);
                % recall fun_l(x,wage,kappa,par)
                emp_loss_x(x_c)  = emp_loss_x(x_c) + ...
                    fun.fun_l(par.A_small(1,i_c)*x_grid(x_c),path_nogrant.w(1),k_grid(k_c),par) ...
                    *dexit*weights_nogrant(k_c,x_c,n_c)*distrib_tran_nogrant.mu(k_c,b_c,x_c,1,n_c);
            end % b_c
        end % k_c
    end % n_c
    exit_rate_x(x_c) = exits_x(x_c)/mass_x(x_c);
end % x_c

% Share of exit due to each x and kappa
exit_share_x = exits_x/sum(exits_x); %dim: (nx,1)

% Save to structure exit_nogrant
exit_nogrant.exit_rate_x  = exit_rate_x;
exit_nogrant.exits_x      = exits_x;
exit_nogrant.exit_share_x = exit_share_x;
exit_nogrant.emp_loss_x   = emp_loss_x;
exit_nogrant.mass_x       = mass_x;
%% Transition with baseline grant, consider only t=1

% State-space: (b,x,kappa,time,impact x grant)

exits_x     = zeros(nx,1); % measure of exits by x (without dividing by mu)
exit_rate_x = zeros(nx,1); % exits/mass of firms equals exit rate
emp_loss_x  = zeros(nx,1); % employment loss due to exit
mass_x      = zeros(nx,1); % mass of incumbent firms by x
for x_c = 1:nx
    for n_c = 1:nn
        [i_c,~] = ind2sub([ni,ns],n_c); % i_c = impact indicator; s_c = grant indicator
        for k_c = 1:nk
            for b_c = 1:nb % current debt
                % dim: (b,x,kappa,time,impact x grant)
                dexit = psi + (1-psi)*pol_tran_grant_baseline.pol_exit(k_c,b_c,x_c,1,n_c);
                mass_x(x_c) = mass_x(x_c) + weights_grant_baseline(k_c,x_c,n_c)*distrib_tran_grant_baseline.mu(k_c,b_c,x_c,1,n_c);

                exits_x(x_c) = exits_x(x_c)+dexit*weights_grant_baseline(k_c,x_c,n_c)*distrib_tran_grant_baseline.mu(k_c,b_c,x_c,1,n_c);
                % recall fun_l(x,wage,kappa,par)
                emp_loss_x(x_c)  = emp_loss_x(x_c) + ...
                    fun.fun_l(par.A_small(1,i_c)*x_grid(x_c),path_grant_baseline.w(1),k_grid(k_c),par) ...
                    *dexit*weights_grant_baseline(k_c,x_c,n_c)*distrib_tran_grant_baseline.mu(k_c,b_c,x_c,1,n_c);
            end % b_c
        end % k_c
    end % n_c
    exit_rate_x(x_c) = exits_x(x_c)/mass_x(x_c);
end % x_c

% Share of exit due to each x and kappa
exit_share_x = exits_x/sum(exits_x); %dim: (nx,1)

% Save to structure exit_nogrant
exit_grant_baseline.exit_rate_x  = exit_rate_x;
exit_grant_baseline.exits_x      = exits_x;
exit_grant_baseline.exit_share_x = exit_share_x;
exit_grant_baseline.emp_loss_x   = emp_loss_x;
exit_grant_baseline.mass_x       = mass_x;


end %END FUNCTION <plot_micro>

