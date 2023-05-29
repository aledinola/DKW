function [agg] = fun_aggregates(par,sol,distribS,phi_dist,prices)
% Usage: [agg] = fun_aggregates(par,sol,distribS,phi_dist,prices);
% Purpose: compute aggregate variables (C,N,etc.) for the steady-state

%------------------------- input checks ----------------------------------%
if isstruct(par)==0
    error("Input <par> in fun_aggregates must be a structure!")
end
if isstruct(sol)==0
    error("Input <sol> in fun_aggregates must be a structure!")
end
% validateattributes(distribS.mu, {'double'}, {'finite', 'nonnan', 'nonempty', 'real',...
%     '>=', 0,'size', [par.nb,par.nx]})
% validateattributes(distribS.mu_active, {'double'}, {'finite', 'nonnan', 'nonempty',...
%     'real', '>=', 0,'size', [par.nb,par.nx]})
% validateattributes(phi_dist, {'double'}, {'finite', 'nonnan', 'nonempty',...
%     'real', '>=', 0,'size', [par.nb,par.nx]})
if isstruct(prices)==0
    error("Input <prices> in fun_aggregates must be a structure!")
end
%-------------------------------------------------------------------------%

% Unpack SOL:
%Note: (1-psi)(1-d^l)=1-(psi+(1-psi)d^l)
pol_exit  = sol.pol_exit;
pol_entry = sol.pol_entry;
pol_kp_ind = sol.pol_kp_ind;
pol_kp     = sol.pol_kp;
% Unpack distribS
mu        = distribS.mu;
mu_active = distribS.mu_active;
entry_vec = distribS.entry_vec;
% Unpack PAR:
%theta    = par.theta;
delta_k  = par.delta_k;
x_grid   = par.x_grid;
k_grid   = par.k_grid;
mass     = par.mass;
psi      = par.psi;
%fixcost  = par.fixcost;
%cost_e   = par.cost_e;
% Unpack PRICES:
wage     = prices.wage;
KL_ratio = prices.KL_ratio;

C_agg = fun.C_foc_labor(wage,par); %Aggregate consumption, from FOC

% Marginal distribution of kappa, mu(k,b,x)
mu_kappa = sum(mu_active,[2,3]);

% Mass of small firms
mass_small = sum(mu_active(:));
% Capital in small firms
K_small = dot(k_grid,mu_kappa);

% Measure of active firms and entrants
Mactive = sum(mu_active(:));
Mentr   = sum(entry_vec(:));

% We include exogenous esit
exit_all = psi+(1-psi)*pol_exit; %dim: (nk,nb,nx)

% Compute aggregate values for entry, output (small firms), liquidation
% (small firms) and employment (small firms)
A = 1; % productivity shock equals 1 in the steady-state
[entry_cost,output_small,liq,L_small,cost_adj,exit,entry] = sub_aggregates_onestep(mu,mu_active,...
     pol_kp_ind,pol_entry,exit_all,phi_dist,wage,mass,A,par);

% Normalize by measure of incumbent firms to get the exit rate:
exit_rate = exit/sum(mu(:));
% Normalize by measure of active firms to get the entry rate:
entry_rate = entry/sum(mu_active(:));


% Compute capital adjustment by type of adjustment

[capadj] = compute_cap_adj(pol_kp,mu,mu_active,phi_dist, ...
    pol_entry,pol_exit,mass,k_grid,delta_k,psi);

% Compute exits excluding microbusinesses
exits_temp = 0;
mass_temp = 0;
for x_c = 1:par.nx
    % We do not multiple with A_small here.
    x_val = x_grid(x_c);
    for k_c = 1:par.nk
        kappa = k_grid(k_c);
        l_opt        = fun.fun_l(x_val,wage,kappa,par);
        if l_opt > par.emp_min
            exits_temp = exits_temp + sum(exit_all(k_c,:,x_c).*mu(k_c,:,x_c),'all');
            mass_temp = mass_temp + sum(mu(k_c,:,x_c),'all');
        end
    end
end
exit_emp = exits_temp;
exit_rate_emp = exits_temp/mass_temp;

% Left-hand side of market clearing eq. on page 49:
LHS = C_agg-output_small+cost_adj+entry_cost-liq;
aux = fun.prod_corp(KL_ratio,1/KL_ratio,par)-delta_k;

% Capital in corporate sector:
K_corp = LHS/aux;
% Employment in corporate sector:
L_corp = K_corp/KL_ratio;
% Output in corporate sector:
Y_corp = fun.prod_corp(KL_ratio,L_corp,par);%A*K_corp^alpha*L_corp^(1-alpha);

% Aggregate employment: L in corporates plus L in small firms
L_agg = L_corp+L_small;

% Y_small is output in small firms, including liq and entry and adj costs:
Y_small = output_small-cost_adj+liq-entry_cost;

% Aggregate output
Y_agg = Y_corp+output_small;

% Capital owned by households (rented to Corporate and small firms)
K_agg = K_corp + K_small;

% Household capital investment in steady-state
InvK_corp = delta_k * K_corp;


% All investment in the economy (\tilde{I} in paper)
InvK = InvK_corp + entry_cost - liq + cost_adj;


% Pack aggregate variables into a structure
agg = v2struct(C_agg,K_agg,K_corp,K_small,mass_small,...
    L_agg,L_corp,L_small,Y_agg,Y_corp,Y_small,output_small,liq,entry,entry_rate,...
    entry_cost,Mactive,Mentr,InvK,InvK_corp,exit_rate,exit,cost_adj,capadj,...
    exit_emp,exit_rate_emp);

end %END FUNCTION <fun_aggregates>

