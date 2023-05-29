function [irf_impact] = fun_aggregates_tran_impact(par,pol_tran,distrib_tran,path,agg_ss)

% #VC# V41

% Purpose: compute aggregate variables (L,K,Y) along transition by impact
% status
% Note: I use _small and _corp to denote var for non-corp and corp sectors
% respectively
%------------------------- input checks ----------------------------------%
if isstruct(par)==0 
    error("Input <par> in fun_aggregates_tran must be a structure!")
end
if isstruct(pol_tran)==0 
    error("Input <pol_tran> in fun_aggregates_tran must be a structure!")
end
if isstruct(distrib_tran)==0 
    error("Input <distrib_tran> in fun_aggregates_tran must be a structure!")
end
if isstruct(path)==0
    error("Input <path> in fun_aggregates_tran must be a structure!")
end
if isstruct(agg_ss)==0
    error("Input <agg_ss> in fun_aggregates_tran must be a structure!")
end

validateattributes(distrib_tran.mu, {'double'}, {'finite', 'nonnan', 'nonempty', 'real',...
    '>=', 0,'size', [par.nb,par.nx,par.nk,par.T+1,par.nn]})
validateattributes(distrib_tran.mu_active, {'double'}, {'finite', 'nonnan', 'nonempty',...
    'real', '>=', 0,'size', [par.nb,par.nx,par.nk,par.T+1,par.nn]})
%-------------------------------------------------------------------------%

% Unpack policies:
pol_exit  = pol_tran.pol_exit;  %dim: (nb,nx,nk,T+1,nn)

% Unpack distributions
mu        = distrib_tran.mu;          %dim: (nb,nx,nk,T+1,nn)
mu_active = distrib_tran.mu_active;   %dim: (nb,nx,nk,T+1,nn)

%mu_ss = distrib_ss.mu; %dim: (nb,nx,nk) 
%mu    = repmat(mu_ss,1,1,1,par.T+1,par.ns);
%mu_active_ss = distrib_ss.mu_active; %dim: (nb,nx,nk) 
%mu_active    = repmat(mu_active_ss,1,1,1,par.T+1,par.ns);

% Unpack PATH for prices
wage_path = path.w;

% Unpack PAR:
T        = par.T;
psi      = par.psi;
fixp     = par.fixp;
k_grid   = par.k_grid;  %dim: (nk,1)
fixcost  = par.fixcost; %dim: (nk,1)
x_grid   = par.x_grid;
nx       = par.nx;
%nb       = par.nb;
nk       = par.nk;
ni       = par.ni;
ns       = par.ns;
A_small  = par.A_small; % dim (T+1,ni)
eta_i    = par.eta_i;
%-------------------------------------------------------------------------%

if par.verbose>=1
    fprintf(' \n')
    disp('Calculate small-firm aggregates on the transition path...')
end

%% Calculate small-firm aggregates on the transition path

% Define weigths 
% - n_c = 1: impact, grant
% - n_c = 2: no impact, grant
% - n_c = 3: impact, no grant
% - n_c = 4: no impact, no grant
weights = par.weights;

% 1. Precompute mass of small firms and labor demand
% Recall: mu in the draft ==> mu_active in the code
%         mu^0 in the draft ==> mu in the code
mass_small =  zeros(T+1,ni);   % mass of small firms
mass_small_0 =  zeros(T+1,ni); % mass of incumbent small firms
for t = 1:T+1
    for i_c = 1:ni % i_c = impact indicator;
        for s_c = 1:ns % s_c = grant indicator;
            n_c = sub2ind([ni,ns],i_c,s_c);
            for k_c = 1:nk
                for x_c = 1:nx
                    % dim: (b,x,k,time,impact x grant)
                    mass_small(t,i_c) = mass_small(t,i_c) + weights(x_c,k_c,n_c)*sum(mu_active(:,x_c,k_c,t,n_c));
                    mass_small_0(t,i_c) = mass_small_0(t,i_c)+ weights(x_c,k_c,n_c)*sum(mu(:,x_c,k_c,t,n_c));
                end
            end
        end
    end
end

% Precompute labor demand function for speed reasons
labor_demand = zeros(nx,nk,T+1,ni);
for i_c = 1:ni % impact 
    for t = 1:T+1
        wage   = wage_path(t);
        for k_c = 1:nk
            kappa = k_grid(k_c);
            for x_c = 1:nx
                x_val = A_small(t,i_c)*x_grid(x_c);
                labor_demand(x_c,k_c,t,i_c) = fun.fun_l(x_val,wage,kappa,par);
            end
        end
    end
end

% 2. Compute employment and output by impacted status

exit_rate_vec = zeros(T+1,ni);
output_small  = zeros(T+1,ni);
L_small       = zeros(T+1,ni);
K_small_owned = zeros(T+1,ni);

for t = 1:T+1
    for i_c = 1:ni % i_c = impact indicator;
        output   = 0;
        capowned = 0;
        exit_rate = 0;
        empl     = 0;
        for s_c = 1:ns %  s_c = grant indicator
            n_c = sub2ind([ni,ns],i_c,s_c);
            for k_c = 1:nk
                kappa = k_grid(k_c); c = fixcost(k_c);
                for x_c = 1:nx
                    x_val = A_small(t,i_c)*x_grid(x_c);
                    empl   = empl + weights(x_c,k_c,n_c)*labor_demand(x_c,k_c,t,i_c)*sum(mu_active(:,x_c,k_c,t,n_c));
                    output = output + weights(x_c,k_c,n_c)*fun.prod_small(x_val,kappa,labor_demand(x_c,k_c,t,i_c),c,par)*sum(mu_active(:,x_c,k_c,t,n_c));
                    capowned = capowned + weights(x_c,k_c,n_c)*(1-fixp)*kappa*sum(mu_active(:,x_c,k_c,t,n_c));
                    psi_tilde = psi+(1-psi)*pol_exit(:,x_c,k_c,t,n_c); % dim (nb,1)
                    exit_rate = exit_rate+weights(x_c,k_c,n_c)*dot(psi_tilde,mu(:,x_c,k_c,t,n_c));
                end %x
            end %k
        end %s_c
        exit_rate_vec(t,i_c) = exit_rate/mass_small_0(t,i_c); % Exit rate
        output_small(t,i_c) = output;           % output small firms
        L_small(t,i_c)      = empl;             % employment small firms
        K_small_owned(t,i_c) = capowned;        % capital owned by small firms
    end %i_c
end %t

%% Calculate irf for mass_small, exit_rate, output_small, L_small, and K_small_owned

impact_weight = [eta_i,1-eta_i];
K_small_owned_ss = agg_ss.K_all - agg_ss.K_agg;


irf_impact = struct();
irf_impact.output_small     = zeros(T+1,ni);
irf_impact.L_small          = zeros(T+1,ni);
irf_impact.K_small_owned    = zeros(T+1,ni);
irf_impact.exit_rate        = zeros(T+1,ni);
irf_impact.mass_small       = zeros(T+1,ni);

for i_c = 1:ni
    irf_impact.mass_small(:,i_c) = (mass_small(:,i_c)-impact_weight(i_c)*agg_ss.mass_small) ...
        /(impact_weight(i_c)*agg_ss.mass_small); % mass of small firms
    irf_impact.output_small(:,i_c) = (output_small(:,i_c)-impact_weight(i_c)*agg_ss.output_small) ...
        /(impact_weight(i_c)*agg_ss.output_small); % output small firms
    irf_impact.L_small(:,i_c) = (L_small(:,i_c)-impact_weight(i_c)*agg_ss.L_small) ...
        /(impact_weight(i_c)*agg_ss.L_small); % employment small firms
    irf_impact.K_small_owned(:,i_c) = (K_small_owned(:,i_c)-impact_weight(i_c)*K_small_owned_ss) ...
        /(impact_weight(i_c)*K_small_owned_ss); % capital small firms
    irf_impact.exit_rate(:,i_c) = (exit_rate_vec(:,i_c)-agg_ss.exit_rate) ...
        /agg_ss.exit_rate; % exit rate small firms
end


end %END FUNCTION <fun_aggregates_tran_impact>

