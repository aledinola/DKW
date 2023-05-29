function [agg_tran] = fun_aggregates_tran(par,pol_tran,distrib_tran,path,agg_ss,distrib_ss,prices_ss)

% #VC# V56

% Purpose: compute aggregate variables (L,K,Y) along transition
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
    '>=', 0,'size', [par.nk,par.nb,par.nx,par.T+1,par.nn]})
validateattributes(distrib_tran.mu_active, {'double'}, {'finite', 'nonnan', 'nonempty',...
    'real', '>=', 0,'size', [par.nk,par.nb,par.nx,par.T+1,par.nn]})
%-------------------------------------------------------------------------%

% Unpack policies:
pol_exit  = pol_tran.pol_exit;  %dim: (nk,nb,nx,T+1,nn)
pol_entry = pol_tran.pol_entry; %dim: (nk,nb,nx,T+1,nn)
phi_dist  = pol_tran.phi_dist; %dim: (nk,nb,nx,T+1,nn)
pol_kp_ind    = pol_tran.pol_kp_ind;
pol_kp    = pol_tran.pol_kp;
% Unpack distributions
mu        = distrib_tran.mu;          %dim: (nk,nb,nx,T+1,nn)
mu_active = distrib_tran.mu_active;   %dim: (nk,nb,nx,T+1,nn)

% Unpack PATH for prices
wage_path = path.w;
KL_ratio  = path.KL_ratio;
C_path    = path.C;

% Unpack PAR:
T        = par.T;
T_grant  = par.T_grant;
%theta    = par.theta;
k_grid   = par.k_grid;  %dim: (nk,1)
%fixcost  = par.fixcost; %dim: (nk,1)
delta_k  = par.delta_k;
x_grid   = par.x_grid;
mass_vec     = par.mass_vec;
nx       = par.nx;
nb       = par.nb;
nk       = par.nk;
ni       = par.ni;
ns       = par.ns;
nn       = par.nn;
A_small  = par.A_small; % dim (T+1,ni)
A_corp   = par.A_corp;
%cost_e   = par.cost_e;
psi      = par.psi;
%eta      = par.eta;  % dim (nx,nk), share of firms with grant
%eta_i    = par.eta_i; % scalar, share of impacted firms
%-------------------------------------------------------------------------%

if par.verbose>=1
    fprintf(' \n')
    disp('Aggregates transition..')
end

% Define weigths 
% - n_c = 1: impact, grant
% - n_c = 2: no impact, grant
% - n_c = 3: impact, no grant
% - n_c = 4: no impact, no grant
weights = par.weights; % dim:(nk,nx,nn)

% 1. Compute demand for rental capital by small firms
% Recall: mu in the draft ==> mu_active in the code
%         mu^0 in the draft ==> mu in the code
K_small = zeros(T+1,1);
mass_small =  zeros(T+1,1);   % mass of small firms
mass_small_0 =  zeros(T+1,1); % mass of incumbent small firms

weights_arr = repmat(weights,1,1,1,nb);
weights_arr = permute(weights_arr,[1 4 2 3]);
kappa_arr = repmat(k_grid(:),1,nb,nx,nn);
mu_active_p = permute(mu_active,[1 2 3 5 4]);
mu_p = permute(mu,[1 2 3 5 4]);


for t = 1:T+1
    K_small(t) = sum(kappa_arr.*weights_arr.*mu_active_p(:,:,:,:,t),'all');
    mass_small(t) = sum(weights_arr.*mu_active_p(:,:,:,:,t),'all');
    mass_small_0(t) = sum(weights_arr.*mu_p(:,:,:,:,t),'all');
end %t


% 2. Compute some aggregate variables in the small firms sector

% Compute capital adjustment by type of adjustment


output_small  = zeros(T+1,1);
Y_small       = zeros(T+1,1); %equals output_small-cost_adj+liq-entry;
entry_vec     = zeros(T+1,1);
entry_rate_vec     = zeros(T+1,1);
entry_cost_vec = zeros(T+1,1);
liq_vec       = zeros(T+1,1);
exit_vec      = zeros(T+1,1);
exit_rate_vec = zeros(T+1,1);
L_small       = zeros(T+1,1);
cost_adj      = zeros(T+1,1);
capadj        = zeros(T+1,4); % four types of capital adjustment, see function compute_cap_adj
exit_emp_vec  = zeros(T+1,1); % exits excluding microbusinesses
exit_rate_emp_vec = zeros(T+1,1); % exit rate excluding microbusinesses
wage_t0       = prices_ss.wage; % wage in the previous period. For t=1, previous period is ss.
A_small_t0    = ones(2,1);
for t = 1:T+1
    % mass of potential entrants at t
    mass = mass_vec(t);
    % reset temporary variables
    output     = zeros(nn,1);
    c_adj      = zeros(nn,1);
    entry      = zeros(nn,1);
    entry_cost = zeros(nn,1);
    exit_rate  = zeros(nn,1);
    liq        = zeros(nn,1);
    empl       = zeros(nn,1);
    exits_emp  = zeros(nn,1);
    mass_emp   = zeros(nn,1);
    capadj_vec = zeros(4,nn);
    
    for n_c = 1:nn
        [i_c,~] = ind2sub([ni,ns],n_c); % i_c = impact indicator; s_c = grant indicator
        
        weights_n  = weights(:,:,n_c);
        weights_nn = permute(repmat(weights_n,1,1,nb),[1,3,2]);
        mu_temp = weights_nn.*mu(:,:,:,t,n_c);
        mu_active_temp = weights_nn.*mu_active(:,:,:,t,n_c);
        phi_dist_temp = weights_nn.*phi_dist(:,:,:,t,n_c);
        % exit_all_n: overall exit rate including exog. and end. exits
        exit_all_n = psi+(1-psi)*pol_exit(:,:,:,t,n_c);  
        % small firm agg variables
        [entry_cost(n_c),output(n_c),liq(n_c),empl(n_c),c_adj(n_c),exit_rate(n_c),entry(n_c)] = sub_aggregates_onestep(mu_temp,...
            mu_active_temp,pol_kp_ind(:,:,:,t,n_c),pol_entry(:,:,:,t,n_c),...
            exit_all_n,phi_dist_temp,wage_path(t),mass,...
            A_small(t,i_c),par);
        % capital adjustment variables
        [capadj_vec(:,n_c)] = compute_cap_adj(pol_kp(:,:,:,t,n_c),mu_temp,mu_active_temp,phi_dist_temp, ...
            pol_entry(:,:,:,t,n_c),pol_exit(:,:,:,t,n_c),mass,k_grid,delta_k,psi);
        % exit rate excluding microbusinesses
        exits_temp = 0;
        mass_temp = 0;
        for x_c = 1:nx
            % A_small and wage from the previous period because we want to know the
            % employment size of the firm in the previous period.
            
            %x_val = A_small_t0(i_c)*x_grid(x_c);
            x_val = x_grid(x_c);
            for k_c = 1:nk
                kappa = k_grid(k_c);
                %l_opt        = fun.fun_l(x_val,wage_t0,kappa,par);
                l_opt        = fun.fun_l(x_val,prices_ss.wage,kappa,par);

                if l_opt > par.emp_min
                    exits_temp = exits_temp + sum(exit_all_n(k_c,:,x_c).*mu_temp(k_c,:,x_c),'all');
                    mass_temp = mass_temp + sum(mu_temp(k_c,:,x_c),'all');
                end
            end
        end
        exits_emp(n_c) = exits_temp;
        mass_emp(n_c) = mass_temp;
    end %n_c
    entry_cost_vec(t) = sum(entry_cost); % cost of entry
    entry_vec(t)      = sum(entry);      % measure of entrants
    entry_rate_vec(t) = sum(entry)/mass_small(t);
    liq_vec(t)      = sum(liq);
    cost_adj(t)     = sum(c_adj);
    exit_vec(t)      = sum(exit_rate);   % measure of exiting firms
    exit_rate_vec(t) = sum(exit_rate)/mass_small_0(t);
    output_small(t) = sum(output);           % output small firms
    L_small(t)      = sum(empl);             % employment small firms
    % Y_small is small firm output minus depreciation plus exit minus entry
    Y_small(t)      = output_small(t)-cost_adj(t)+liq_vec(t)-entry_cost_vec(t);  
    capadj(t,:)     = sum(capadj_vec,2)';
    % exits excluding microbusinesses
    exit_emp_vec(t) = sum(exits_emp);
    exit_rate_emp_vec(t) = sum(exits_emp)/sum(mass_emp);
    % update A_small and wage in the previous period
    wage_t0 = wage_path(t);
    A_small_t0 = A_small(t,:);
end %t

% 3. Compute labor and output in the corporate sector and investment in K
L_corp = zeros(T+1,1);
Y_corp = zeros(T+1,1);
K_corp = zeros(T+1,1);
%K_agg  = zeros(T+1,1); %equals K_corp + K_small
InvK_corp   = zeros(T+1,1); % Corporate investment

% Capital is predetermined
%K_agg(1)  = agg_ss.K_agg;
K_corp(1) = agg_ss.K_corp; 

for t = 1:T+1
    % two ways to compute K_corp
    % (A) use resource constraint:
    %LHS = C_path(t)+entry-output-liq+delta_k*K_small(t);
    %aux = KL_ratio(t)^(alpha-1)-delta_k;
    %K_corp1(t) = LHS/aux;
    % (B) use definition of aggregate capital:
    %K_corp(t) = K_agg(t)-K_small(t);
    L_corp(t) = K_corp(t)/KL_ratio(t);
    if K_corp(t)>0
        Y_corp(t) = A_corp(t)*fun.prod_corp(KL_ratio(t),L_corp(t),par);
    else
        Y_corp(t) = 0;
    end

    InvK_corp(t)   = Y_small(t)+Y_corp(t)-C_path(t);
    if (t<T+1)
        K_corp(t+1) = (1-delta_k)*K_corp(t)+InvK_corp(t);
    end
end
K_agg = K_corp + K_small;
% if any(K_corp<0)
%     warning("K_corp<0 for some t")
%     keyboard
% end

Y_agg = Y_corp+output_small;
L_agg = L_corp+L_small;

% Economy wide capital and investment
InvK = InvK_corp + entry_cost_vec - liq_vec + cost_adj;

% Cost of grant policy

%pol_tran.grant_vec dim: nx,nk,T,ns
tot_grant = 0;
for i_c = 1:ni
    n_c = sub2ind([ni,ns],i_c,1);
    for t = 1:T_grant
        for k_c = 1:nk
            for x_c = 1:nx
                % all exits
                exit_all_b = psi+(1-psi)*squeeze(pol_exit(k_c,:,x_c,t,n_c));
                tot_grant = tot_grant + pol_tran.grant_vec(k_c,x_c,t,1)*weights(k_c,x_c,n_c) ...
                    *sum(dot(1-exit_all_b,squeeze(mu(k_c,:,x_c,t,n_c))));
            end
        end
    end
end

agg_tran = v2struct(Y_agg,Y_corp,Y_small,output_small,...
    K_small,K_corp,K_agg,InvK,InvK_corp,...
    L_agg,L_corp,L_small,mass_small,mass_small_0,...
    entry_vec,entry_cost_vec,entry_rate_vec,liq_vec,exit_rate_vec,exit_vec,cost_adj,tot_grant,capadj,...
    exit_emp_vec,exit_rate_emp_vec);


end %END FUNCTION <fun_aggregates_tran>

