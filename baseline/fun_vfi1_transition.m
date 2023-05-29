function [polS] = fun_vfi1_transition(par,path,sol_ss,b_grid_ss)
% Purpose: solve firm's value funct iter during the transition
% Usage: [polS] = fun_vfi1_transition(par,path,sol_ss,b_grid_ss)

%%%%% Overview:
%   This function implements steps 2a and 2b outlined in fun_transition.m:
    %   2.a. Calculate profit profit_mat(1:xn,1:T+1)
    %   2.b. Calculate value functions by backward induction val(1:nb,1:nx,1:T+1),val0(1:nb,1:nx,2:T+1);
    %       Also calculate associated policy functions pol_entry,pol_exit,
    %       pol_debt for t = 1,...T+1
%%%%% Value and policy functions in the steady-state:
% val       :: v(b,x)
% pol_debt  :: b'(b,x)
% pol_exit  :: d^l(b,x)
% pol_entry :: d^e(b,x)

%------------------------- input checks ----------------------------------%
% if isstruct(path)==0
%     error("Input <path> in VFI must be a structure!")
% end
% if isstruct(sol_ss)==0
%     error("Input <sol_ss> in VFI must be a structure!")
% end
% if isstruct(par)==0
%     error("Input <par> in VFI must be a structure!")
% end
%-------------------------------------------------------------------------%

% Unpack some flags
verbose   = par.verbose;

%Unpack some parameters:
T       = par.T;
Tmax    = par.Tmax;
theta   = par.theta;
x_grid  = par.x_grid;
pi_x    = par.pi_x;
nx      = par.nx;
nb      = par.nb;
ns      = par.ns;
ni      = par.ni; % impacted and unmipacted
nn      = par.nn; % = ns*ni
nk      = par.nk;
delta   = par.delta_k; %depreciation rate
lambda_vec  = par.lambda_vec;
psi     = par.psi;
cost_e  = par.cost_e;
fixcost = par.fixcost; % dim nk
k_grid  = par.k_grid; % dim nk
A_small = par.A_small; % dim (T+1,2)
%eta_i   = par.eta_i;    % fraction of impacted
Xp      = par.Xp;      % scalar
T_grant = par.T_grant; % number of perios grant is given
%x0_prob = par.x0_prob; % distrib of x for entrants


% Unpack last period values from the steady state:
V1_final      = sol_ss.V1;   % This is V(k,x) for unc firms, dim:(nk,nx)
B_hat_final   = sol_ss.B_hat;  % dim: (nk,nx)
x_tilde_val_final = sol_ss.x_tilde_val; % dim: (nk,1)
val_final     = sol_ss.val; %dim: (nb,nx,nk)
wage_ss       = sol_ss.wage; 
phi_dist_ss   = sol_ss.phi_dist; %dim: (nb,nx,nk)
val0_unc_final = sol_ss.val0_unc;
% Unpack transition paths for prices
q_path      = path.q;
wage_path   = path.w;

assert(isequal(size(A_small),[T+1,ni]),'A_small has wrong dimensions!')

nn_vec             = zeros(T+1,1);
for t = 1:T+1
   if t<=Tmax
       nn_vec(t)=nn;
   else
       nn_vec(t)=1;
   end
end

%%   2.a. Calculate profit profit_mat(1:nk,1:xn,1:T+1,1:2)
% Note: the grant does not change optimal labor demand

grant_vec  = zeros(nk,nx,T+1,ns); %dim: (nk,nx,T+1,ns)
for x_c = 1:nx
    for k_c = 1:nk
        kappa = k_grid(k_c);
        for t = 1:T_grant
            % Grant amount is based on labor demand *if there were NO TFP
            % shock (so, doesn't depend on impactedness)
            x_val = x_grid(x_c);
            grant_vec(k_c,x_c,t,1) = Xp*wage_ss*fun.fun_l(x_val,wage_ss,kappa,par)/T_grant;
        end
    end
end

profit_mat = zeros(nk,nx,T+1,nn);
for n_c = 1:nn % 1 = impact, grant,  2 = no impact, grant, 3 = impact, no grant, 4 =  no impact, no grant,
    [i_c,s_c] = ind2sub([ni,ns],n_c); % i_c = impact indicator; s_c = grant indicator
    for t = 1:T+1
        for x_c = 1:nx
            x_val = x_grid(x_c)*A_small(t,i_c);
            for k_c = 1:nk
                kappa = k_grid(k_c);
                c     = fixcost(k_c);               
                profit_mat(k_c,x_c,t,n_c) = fun.fun_profit(x_val,kappa,c,wage_path(t),par);
                % Add grant
                profit_mat(k_c,x_c,t,n_c) = profit_mat(k_c,x_c,t,n_c)+grant_vec(k_c,x_c,t,s_c);
            end
        end
    end
end

%%   2.b. Calculate value functions by backward induction 


%% VFI for unconstrained firms by backward induction
if verbose>=1; disp('VFI for unconstrained firms...'); end

V1         = zeros(nk,nx,T+1,nn);
Vnext      = repmat(V1_final,1,1,nn); % Terminal condition at T+1 is the steady-state

for t=T+1:-1:1
    for n_c = 1:nn_vec(t)
        V1(:,:,t,n_c) = sub_V1_onestep(Vnext(:,:,n_c),pi_x,profit_mat(:,:,t,n_c), ...
            k_grid,q_path(t),theta,delta,psi,0);
        % Update
        Vnext(:,:,n_c) = V1(:,:,t,n_c);
        if t>Tmax
            for nt_c = 2:nn
                V1(:,:,t,nt_c)  = V1(:,:,t,1);
                Vnext(:,:,nt_c) = Vnext(:,:,1);
            end
        end
    end % t
end % n_c

if verbose>=2
    fprintf(" \n")
end

validateattributes(V1, {'double'}, {'finite', 'nonnan', 'nonempty', 'real','size', [nk,nx,T+1,nn]})
if verbose>=1; disp("VFI for unconstrained firms done!"); end

%% Find the productivity cutoff x_tilde

x_tilde     = zeros(nk,T+1,nn);
x_tilde_val = zeros(nk,T+1,nn);
for n_c = 1:nn
    for t=T+1:-1:1
        for k_c = 1:nk
            k_val    = k_grid(k_c);
            viable_x = find(squeeze(V1(k_c,:,t,n_c))>=theta*(1-delta)*k_val);
            if ~isempty(viable_x)
                x_tilde(k_c,t,n_c)     = min(viable_x);
                x_tilde_val(k_c,t,n_c) = x_grid(x_tilde(k_c,t,n_c));
            else
                fprintf('There are no viable x for k_c=%d, k_val=%f, t=%d, n_c=%d \n',...
                    k_c,k_val,t,n_c)
                x_tilde(k_c,t,n_c)     = nx+1;
                x_tilde_val(k_c,t,n_c) = x_grid(nx)+1;
            end %if
        end % k_c
    end %t
end %n_c


%% Compute unconstrained SS investment policy

Vnext = repmat(V1_final,1,1,nn); % Terminal condition at T+1 is the steady-state
pol_kp_unc = zeros(nk,nx,T+1,nn);
for n_c = 1:nn
    for t=T+1:-1:1
        [pol_kp_unc(:,:,t,n_c)] = sub_investment_onestep(Vnext(:,:,n_c),...
            k_grid,pi_x,theta,delta,q_path(t),psi);
        Vnext(:,:,n_c) = V1(:,:,t,n_c);
    end
end

%% Compute B_hat(k,x) by backward induction

if verbose>=1; disp("Fixed point B_hat..."); end
% Allocate B_hat(k,x,t,n) and b'(k,x,t,n) unconstrained
B_hat      = ones(nk,nx,T+1,nn);
pol_bp_unc = zeros(nk,nx,T+1,nn);
B_hat_next = repmat(B_hat_final,1,1,nn); 
x_tilde_val_next = repmat(x_tilde_val_final,1,nn); 
for n_c = 1:nn
    for t=T+1:-1:1
        [B_hat(:,:,t,n_c),pol_bp_unc(:,:,t,n_c)] = sub_Bhat_onestep(B_hat_next(:,:,n_c),...
            pol_kp_unc(:,:,t,n_c),profit_mat(:,:,t,n_c),x_tilde_val_next(:,n_c),k_grid,...
            x_grid,q_path(t),theta,delta,lambda_vec(t));
        B_hat_next(:,:,n_c)     = B_hat(:,:,t,n_c);
        x_tilde_val_next(:,n_c) = x_tilde_val(:,t,n_c);
    end %t
end %n_c

%% Define flexible grid for b, dependent on k

b_tilde = zeros(nk,T+1,nn); % b_tilde = min(B_hat) for each k
b_grid  = zeros(nb,nk,T+1,nn);
for n_c = 1:nn
    for t=1:T+1
        for k_c=1:nk
            k_val = k_grid(k_c);
            viable_x = x_grid>=x_tilde_val(k_c,t,n_c);
            if sum(viable_x)==0
                fprintf('k_c = %d, k_val = %f \n',k_c,k_val)
                warning('There are NO viable x''s ')
                b_tilde(k_c,t,n_c) = lambda_vec(t)*k_val;
                b_grid(:,k_c,t,n_c)= lambda_vec(t)*k_val;
            else
                b_tilde(k_c,t,n_c)  = min(squeeze(B_hat(k_c,viable_x,t,n_c)));
                b_grid(:,k_c,t,n_c) = linspace(min(b_tilde(k_c,t,n_c),lambda_vec(t)*k_val),lambda_vec(t)*k_val,nb)';
            end
        end
    end
end
b_grid = permute(b_grid,[2,1,3,4]);

%% Precompute upper bound for next-period capital k'

% Upper bound depends on whether the capital adjustment is upward or
% downward
kp_bar = zeros(nk,nb,nx,T+1,nn);
for n_c = 1:nn
    for t=T+1:-1:1
        kp_bar(:,:,:,t,n_c) = sub_kp_onestep(profit_mat(:,:,t,n_c),b_grid(:,:,t,n_c),...
            k_grid,q_path(t),theta,delta,lambda_vec(t));
    end
end

%% Value function v(b,x)=V(x)-b for unconstrained firms

% Value function unconstrained firms: v(b,x,k) = V(x,k)-b
val_unc  = zeros(nk,nb,nx,T+1,nn); % after exit (i.e. continuing firms)
val0_unc = zeros(nk,nb,nx,T+1,nn); % before exit
for n_c = 1:nn
    for t=1:T+1
        for x_c = 1:nx
            for b_c = 1:nb
                for k_c = 1:nk
                    k_val = k_grid(k_c);
                    b_val = b_grid(k_c,b_c,t,n_c);
                    val_unc(k_c,b_c,x_c,t,n_c)  = V1(k_c,x_c,t,n_c)-b_val;
                    val0_unc(k_c,b_c,x_c,t,n_c) = max(theta*(1-delta)*k_val-b_val,...
                        val_unc(k_c,b_c,x_c,t,n_c));
                end
            end
        end
    end
end

%% Solve problem of constrained firms

val = zeros(nk,nb,nx,T+1,nn); % v_t(k,b,x,t,imp x grant)
pol_kp_ind_con = ones(nk,nb,nx,T+1,nn);
% val_next: terminal condition for value after entry
val_next       = repmat(val_final,1,1,1,nn); % dim:(nk,nb,nx,nn)
val0_unc_next  = repmat(val0_unc_final,1,1,1,nn); % dim:(nk,nb,nx,nn)

% b_tilde_next in the final period is the ss b_tilde
% same for b_grid_next
b_grid_next = repmat(b_grid_ss,1,1,nn); % dim:(nk,nb,nn)
B_hat_next  = repmat(B_hat_final,1,1,nn); %dim:(nk,nn)

if verbose>=1; disp("VFI for constrained firms..."); end
for t=T+1:-1:1
 
    for n_c = 1:nn_vec(t)
        if verbose>=1
            fprintf('t = %d ; n_c = %d \n',t,n_c)
        end
        [val(:,:,:,t,n_c),pol_kp_ind_con(:,:,:,t,n_c)] = sub_vfi_onestep_mex(val_next(:,:,:,n_c),...
            val0_unc_next(:,:,:,n_c),kp_bar(:,:,:,t,n_c),B_hat_next(:,:,n_c),...
            profit_mat(:,:,t,n_c),k_grid,b_grid_next(:,:,n_c),pi_x,theta,q_path(t),...
            delta,psi,int32(0),int32(0));
        if t>Tmax
            for nt_c = 2:nn
                val(:,:,:,t,nt_c) = val(:,:,:,t,1);
                pol_kp_ind_con(:,:,:,t,nt_c) = pol_kp_ind_con(:,:,:,t,1);
            end
        end
    end %n_c
    % Update
    val_next      = squeeze(val(:,:,:,t,:));
    val0_unc_next = squeeze(val0_unc(:,:,:,t,:));
    b_grid_next   = squeeze(b_grid(:,:,t,:));
    B_hat_next    = squeeze(B_hat(:,:,t,:));
end %t
       
if verbose>=1
    fprintf(' \n')
end

validateattributes(val, {'double'}, {'finite', 'nonnan', 'nonempty', 'real','size', [nk,nb,nx,T+1,nn]})

%% Update pol_debt, pol_kp, pol_kp_ind and val
% val is the value of all firms, constrained and unconstrained 
pol_debt   = zeros(nk,nb,nx,T+1,nn);
pol_kp     = zeros(nk,nb,nx,T+1,nn);
pol_kp_ind = ones(nk,nb,nx,T+1,nn);

for n_c = 1:nn
    for t=T+1:-1:1
        [pol_debt(:,:,:,t,n_c),pol_kp(:,:,:,t,n_c),pol_kp_ind(:,:,:,t,n_c),val(:,:,:,t,n_c)] = ...
            fun_pol_update(val(:,:,:,t,n_c),val_unc(:,:,:,t,n_c),pol_bp_unc(:,:,t,n_c),...
            pol_kp_unc(:,:,t,n_c),pol_kp_ind_con(:,:,:,t,n_c),profit_mat(:,:,t,n_c),...
            B_hat(:,:,t,n_c),k_grid,b_grid(:,:,t,n_c),q_path(t),theta,delta);
    end
end

%% Entry and exit policies

% Exit policy d^l(b,x,kappa,time,impact x grant) 
pol_exit  = zeros(nk,nb,nx,T+1,nn);
pol_entry = zeros(nk,nb,nx,T+1,nn);
for t=1:T+1
    for n_c = 1:nn_vec(t)
        % Find entry and exit policy
        [pol_entry(:,:,:,t,n_c),pol_exit(:,:,:,t,n_c)] = ...
            fun_entry_exit(val(:,:,:,t,n_c),profit_mat(:,:,t,n_c),b_grid(:,:,t,n_c),...
            k_grid,theta,delta,cost_e);
        % Interpolate
        [pol_entry(:,:,:,t,n_c),pol_exit(:,:,:,t,n_c)] = ...
            interp_entry_exit(pol_entry(:,:,:,t,n_c),pol_exit(:,:,:,t,n_c), ...
            b_grid(:,:,t,n_c),val(:,:,:,t,n_c),profit_mat(:,:,t,n_c),par);
        % Make copies to other n_c's
        if t>Tmax
            for nt_c = 2:nn
                pol_entry(:,:,:,t,nt_c) = pol_entry(:,:,:,t,1);
                pol_exit(:,:,:,t,nt_c) = pol_exit(:,:,:,t,1);
            end
        end %if t>Tmax  
    end % t
end % n_c


% Entrants do not receive grant, so we impose s=2 (i.e. no grant) on all
% cases.
for n_c = 1:nn
    [i_c,s_c] = ind2sub([ni,ns],n_c); % i_c = impact indicator; s_c = grant indicator
    if s_c == 1
        % find the index nx_c such that i_c is the same and grant indicator
        % is no grant
        nx_c=sub2ind([ni,ns],i_c,ns);
        pol_entry(:,:,:,:,n_c) = pol_entry(:,:,:,:,nx_c);
    end
end %n_c
validateattributes(pol_exit, {'double'}, {'finite', 'nonnan', 'nonempty','>=',0,'<=',1,'real','size', [nk,nb,nx,T+1,nn]})
validateattributes(pol_entry, {'double'}, {'finite', 'nonnan', 'nonempty','real','>=',0,'<=',1,'size', [nk,nb,nx,T+1,nn]})

%% Initial distribution on the transition path

% Note: since entrants do not receive the PPP loan, phi_dist doesn't depend
% on loan or x0
% phi_dist has dim: (nb,nx,nk,T+1,nn)
[phi_dist] = fun_phi_tran(phi_dist_ss,b_grid_ss,b_grid,b_tilde,par);

%% Pack solutions into a structure

polS = v2struct(val,pol_debt,pol_kp_ind,pol_kp,pol_kp_unc,pol_exit,pol_entry,b_tilde,B_hat,b_grid,phi_dist,grant_vec,V1);


end %END FUNCTION <fun_vfi1_transition>
