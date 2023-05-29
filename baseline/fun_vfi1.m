function [sol,b_grid,phi_dist,conv_flag] = fun_vfi1(prices,par)
% Purpose: solve firm's dynamic programming problem, given prices {q,w,R}
% Usage: [sol,b_grid,phi_dist,conv_flag] = fun_vfi(prices,par)
% where the inputs prices and par are structures
%%%%% Overview:
%     (1) given prices (q,w,R), compute labor demand function and current-period 
%         profit \pi(x)
%     (2) Given prices and \pi(x), solve Bellman equation for unconstrained
%         firms (20) obtaining V(x) and cutoff x_tilde, b_tilde
%     (3) Define b_grid = [b_tilde:theta*kappa], V_unc(x,b)=V(x)-b. Obtain V(x,b)
%         for constrained firms using eq. (21)-(22)-(23). Obtain exit rules
%         from eq. (13)-(14).
%%%%% Value and policy functions:
% val       :: v(b,x), after exit
% val0      :: v^0(b,x), before exit
% pol_debt  :: b'(b,x), next-period debt
% pol_exit  :: d^l(b,x), liquidation/exit policy
% pol_entry :: d^e(b,x), entry policy
% Here we also compute b_grid and phi(x,b)
% NOTES
% sub_vfi_onestep_mex is a MEX file. It is saved in the folder with Fortran
% source codes. Compile it by typing the following in Matlab:
% mex -v OPTIMFLAGS="/O3 /Qprec-div- /QxHost /DNDEBUG" -R2018a sub_vfi_onestep_mex.f90 -output sub_vfi_onestep_mex

%------------------------- input checks ----------------------------------%
if isstruct(prices)==0
    error('Input <prices> in fun_vfi1 must be a structure!')
end
if isstruct(par)==0
    error('Input <par> in fun_vfi1 must be a structure!')
end
%-------------------------------------------------------------------------%

% Unpack some flags
verbose   = par.verbose;
max_iter  = par.max_iter;
tolerance = par.tol_vfi;
tol_vfi_u = par.tol_vfi_u;
do_howard = par.do_howard;
n_howard  = par.n_howard;
tol_bhat  = par.tol_bhat;

%Unpack some parameters:
theta   = par.theta;
delta   = par.delta_k;
lambda  = par.lambda;
cost_e  = par.cost_e;
x_grid  = par.x_grid;
pi_x    = par.pi_x;
nx      = par.nx;
nb      = par.nb;
nk      = par.nk;
x0_prob = par.x0_prob;
k_grid  = par.k_grid;
fixcost = par.fixcost; % vector (nk,1)
prob_k  = par.prob_k;
psi     = par.psi;

bk0_vec   = par.bk0_vec;
bk0_prob  = par.bk0_prob;

% Unpack prices:
q      = prices.q;
wage   = prices.wage;

if verbose>=1
    fprintf('--------------------------------------------\n')
    fprintf('VFI \n')
    fprintf('--------------------------------------------\n')
end

% Initialize convergence flag
conv_flag = 0;

%% Compute value of unconstrained firms net of debt, i.e. V(x)

% Precompute static profit on the grid for (x,k)
profit_mat = zeros(nk,nx);
for x_c = 1:nx
    for k_c = 1:nk
        profit_mat(k_c,x_c) = fun.fun_profit(x_grid(x_c),k_grid(k_c),fixcost(k_c),wage,par);
    end
end

V1     = zeros(nk,nx);
ind    = 0;
errter = 10000;
if verbose>=1; disp('VFI for unconstrained firms...'); end
tic
while ind<max_iter && errter > tol_vfi_u

    V2 = sub_V1_onestep(V1,pi_x,profit_mat,k_grid,q,theta,delta,psi,do_howard);
    
    errter = max(abs(V1-V2),[],'all');
    ind    = ind+1; 
    % Update
    V1 = V2;
    
    if verbose>=2
        fprintf('iter = %d, err = %f \n',ind,errter)
    end
end
toc

%% Find the productivity cutoff x_tilde

x_tilde     = zeros(nk,1);
x_tilde_val = zeros(nk,1);
for k_c = 1:nk
    k_val    = k_grid(k_c);
    viable_x = find(V1(k_c,:)>=theta*(1-delta)*k_val);
    if ~isempty(viable_x)
        x_tilde(k_c) = min(viable_x);
        x_tilde_val(k_c) = x_grid(x_tilde(k_c));
    else
        fprintf('There are no viable x for k_c=%d, k_val=%f \n',k_c,k_val)
        conv_flag = -1;
        sol=[];b_grid=[];phi_dist=[];
        return
    end
end

if ~all(isfinite(x_tilde_val))
    warning('x_tilde_val is not finite')
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end

%% Compute unconstrained SS investment policy

[pol_kp_unc] = sub_investment_onestep(V1,k_grid,pi_x,theta,delta,q,psi);

%% Compute B_hat(k,x) as fixed point of (21) and (22)

B_hat  = ones(nk,nx);
ind    = 0;
errter = 100;
if verbose>=1; disp('Fixed point B_hat...'); end
tic
while ind<max_iter && errter > tol_bhat
    
    [B_hat_new,pol_bp_unc] = sub_Bhat_onestep(B_hat,pol_kp_unc,profit_mat,...
        x_tilde_val,k_grid,x_grid,q,theta,delta,lambda);
    
    errter = max(abs(B_hat-B_hat_new),[],'all');
    ind = ind+1; 
    
    % Update
    B_hat = B_hat_new;
    
    if verbose>=2
        fprintf('iter = %d, err = %f \n',ind,errter)
    end
end
toc

if errter>tol_bhat
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end

%% Define flexible grid for b, dependent on k

b_tilde = zeros(nk,1);
b_grid  = zeros(nk,nb);
for k_c=1:nk
    k_val = k_grid(k_c);
    viable_x = x_grid>=x_tilde_val(k_c);
    if sum(viable_x)==0
        fprintf('k_c = %d, k_val = %f \n',k_c,k_val)
        warning('There are NO viable x''s ')
        keyboard
    end
    b_tilde(k_c)  = min(B_hat(k_c,viable_x));
    b_grid(k_c,:) = linspace(min(b_tilde(k_c),lambda*k_val),lambda*k_val,nb);
end

%% Precompute upper bound for next-period capital k'

% Upper bound depends on whether the capital adjustment is upward or
% downward (see eq. 25 draft)
kp_bar = sub_kp_onestep(profit_mat,b_grid,k_grid,q,theta,delta,lambda);

validateattributes(kp_bar, {'double'}, {'finite', 'nonnan', 'nonempty', 'real','size', [nk,nb,nx]})

%% Value function v(b,x)=V(x)-b for unconstrained firms

% Value function unconstrained firms: v(b,x,k) = V(x,k)-b
val_unc  = zeros(nk,nb,nx); % after exit (i.e. continuing firms)
val0_unc = zeros(nk,nb,nx); % before exit

for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            k_val = k_grid(k_c);
            b_val = b_grid(k_c,b_c);
            val_unc(k_c,b_c,x_c)  = V1(k_c,x_c)-b_val;
            val0_unc(k_c,b_c,x_c) = max(theta*(1-delta)*k_val-b_val,val_unc(k_c,b_c,x_c));
        end
    end
end

%% Solve problem of constrained firms using eq. (23)-(25)

% Here "val" is v_c(x,b,k), initial guess for value constrained firm after exit
val     = val_unc; 
iter    = 0;
dist    = 100;
if verbose>=1; disp('VFI for constrained firms...'); end

% Fixed point val --> val_new
while iter<max_iter && dist > tolerance
   
    [val_new,pol_kp_ind_con] = sub_vfi_onestep_mex(val,val0_unc,kp_bar,B_hat,...
        profit_mat,k_grid,b_grid,pi_x,theta,q,delta,psi,int32(do_howard),int32(n_howard));
    
    dist = max(abs(val_new(:)-val(:)));
    iter = iter+1; 
    % Update
    val  = val_new; %0.2*val_new+0.8*val;
    if verbose>=2
        fprintf('iter = %d, err = %.18f \n',iter,dist)
    end
end
if dist>tolerance
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end

% Compute optimal debt policy implied by eq. (25)
[pol_debt,pol_kp,pol_kp_ind,val] = fun_pol_update(val,val_unc,pol_bp_unc,...
    pol_kp_unc,pol_kp_ind_con,profit_mat,B_hat,k_grid,b_grid,q,theta,delta);

%% Compute entry and exit policy 

[pol_entry,pol_exit,pol_exit_forced,pol_exit_vol] = ...
    fun_entry_exit(val,profit_mat,b_grid,k_grid,theta,delta,cost_e);

% Interpolate entry/exit
 [pol_entry,pol_exit] = interp_entry_exit(pol_entry,pol_exit, ...
    b_grid,val,profit_mat,par);

if ~all(isfinite(b_grid))
    warning('b_grid is not finite')
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end
if ~all(isfinite(val))
    warning('val is not finite')
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end
if ~all(isfinite(pol_debt))
    warning('pol_debt is not finite')
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end
if ~all(isfinite(pol_exit))
    warning('pol_exit is not finite')
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end
if ~all(isfinite(pol_entry))
    warning('pol_entry is not finite')
    conv_flag = -1;
    sol=[];b_grid=[];phi_dist=[];
    return
end

%% Initial distribution of (k,b,x) for entrants (\Phi in the draft). 
% Recall that the distrib for k,b,x are independent from each other.
% Entrants draw (k,x) from two non-degenerate distributions: prob_k (Pareto)
% and x0_prob (Log-Normal). 
% b0 computed based on initial debt-to-capital ratio par.bk0

[phi_dist] = gen_phi_dist(bk0_vec,bk0_prob,k_grid,b_grid,x0_prob,prob_k);

%% Pack solutions into a structure
sol  = v2struct(V1,val,pol_kp_ind,pol_kp,pol_kp_unc,pol_debt,pol_exit,pol_exit_forced, ...
    pol_exit_vol,pol_entry,val_unc,val0_unc,profit_mat,x_tilde_val,b_tilde,wage,kp_bar,...
    B_hat);

% Output checks
if (isstruct(sol)==0)
    error('Output <sol> must be a structure')
end

end %END FUNCTION <fun_vfi>
