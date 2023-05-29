function [agg_tran,path,conv_flag,pol_tran,distrib_tran] = fun_transition(par,sol_ss,agg_ss,distrib_ss,prices_ss,b_grid_ss)
% This function computes the transition after an unexpected pandemic shock
% at t=0
% INPUTS:
% par        :: struct with parameters
% sol_ss     :: struct solution steady-state (policy functions and values)
% agg_ss     :: struct with aggregates in the steady-state
% distrib_ss :: struct with distributions in the steady-state
% distrib_ss :: struct with prices in the steady-state

%{ 
             IMPLEMENTATION OF IMPACTED VS UNIMPACTED
for every period on the trantion path, we solve for the value functions and
mu's (and mu_active's) separately for "impacted" and "unimpacted" firms.
In calculating the aggregate measures (e.g. K_small, L_small, Y_small), we
calculated to weighted sum of impacted and unemployed firms using eta and
1-eta. 
Note that the policy functions (entry, exit, debt) in periods t>1 should be
identical for impacted and unimpacted firms, but there would be a
discrepency in the mu's.
%}

%% Unpack Parameters
T           = par.T; % number of transition periods
max_iter_tr = par.max_iter_tr;
tol_tran    = par.tol_tran;
damp        = par.dampening;
verbose     = par.verbose;
disp_tran   = par.disp_tran;

beta        = par.beta;
A_corp      = par.A_corp;
margutil    = par.margutil;
lsupply     = par.lsupply;
delta_k     = par.delta_k;
zeta        = par.zeta;
sigma       = par.sigma;

%% Initialize Paths
path.C        = zeros(T+1,1);
path.KL_ratio = zeros(T+1,1);
path.w        = zeros(T+1,1);
path.q        = zeros(T+1,1);

%%   0. Guess consumption in period 1 (C_path(1))

% Period t=0 in the draft is period t=1 in the code
C_h    = vpa(0.85*agg_ss.C_agg);
C_l    = vpa(0.95*agg_ss.C_agg);
C_init = vpa(0.5*(C_h+C_l));

iter = 0;
err_tran = tol_tran+1;
% Options for fzero
options = optimset('FunValCheck','on','TolX',1e-20);

%% Start transition loop

while abs(err_tran)>tol_tran && iter<=max_iter_tr
    if disp_tran==1
        tic
    end
    
%%   1. Solve for C_path(2:T+1), wage(1:T+1),KL_ratio(1:T+1), q(1:T+1)
   
    path.C(1) = C_init;
    
    for t=1:T+1
            
        % wage(t)
        path.w(t) = lsupply(t)*zeta*(path.C(t))^sigma;

        % KL_ratio(t)
        if t==1
            path.KL_ratio(t) = fun.KL_tran(path.w(t),A_corp(t),par);
        end
        % KL_ratio(t+1)
        if t<T+1
            fun_KL_ratio = @(k_next) fun.dyn_eqn_capital(k_next,path.KL_ratio(t),par,t);
            [path.KL_ratio(t+1),fval,flag] = fzero(fun_KL_ratio,path.KL_ratio(t),options);
            if flag<0
                warning("fzero did not converge")
                fprintf("function value = %f \n", fval)
            end
        end
        % q(t)
        if t<T+1
            %
            path.q(t)=1/(1-delta_k+A_corp(t+1)*fun.marg_prod_capital(path.KL_ratio(t+1),par));
        else
            path.q(t) = prices_ss.q;
        end

        % C(t+1)
        if t<T+1
           path.C(t+1) = path.C(t)*((beta*margutil(t+1))/(path.q(t)*margutil(t)))^(1/sigma);
        end  
       
    end
    
    %% 2. Value function: backward iteration
    %   2.a. Calculate profit profit_mat(1:xn,1:T+1)
    %   2.b. Calculate value functions by backward induction val(1:nb,1:nx,1:T+1),val0(1:nb,1:nx,2:T+1);
    %       Also calcualte associated policy functions pol_entry,pol_exit,
    %       pol_debt for t = 1,...T+1
    if verbose>=1
        disp('--------------------------------------------')
        fprintf('START VFI TRANSITION.. \n')
        tic; 
    end
    [pol_tran] = fun_vfi1_transition(par,path,sol_ss,b_grid_ss);
    if verbose>=1
        time=toc;
        fprintf('Time to do VFI TRANSITION: %8.4f \n', time)
        fprintf(' \n')
    end
    
    %% 3. Distribution: forward iteration

    % Inputs for distribution: steady state mu and mu_active, policy functions
    % for debt, entry and liquidation for all t, computed in VFI
    % Outputs: Sequences of distributions for all t
    if verbose>=1; tic; end
    [distrib_tran] = fun_distrib1_tran(par,pol_tran,distrib_ss,sol_ss);
    if verbose>=1
        time=toc;
        fprintf('Time to do DISTRIBUTION: %8.4f \n', time)
        fprintf(' \n')
    end
    
    %% Aggregates along the transition
    [agg_tran] = fun_aggregates_tran(par,pol_tran,distrib_tran,path,agg_ss,distrib_ss,prices_ss);
    K_agg = agg_tran.K_agg;
    
   %% Update C(1) - find C(1) s.t. agg_tran.K_agg(1) is close to agg_ss.K_agg
   
   if (any(K_agg<0) || any(~isreal(K_agg)))
       % K is too low, decrease C(1)
       C_h = vpa(damp*C_init+(1-damp)*C_h);
       %C_h = C_init;
   else
       if K_agg(T+1)<agg_ss.K_agg
           C_h = vpa(damp*C_init+(1-damp)*C_h);
           %C_h = C_init;
       else
           C_l = vpa(damp*C_init+(1-damp)*C_l);
           %C_l = C_init;
       end
   end
   %% Update
   err_tran_vec = K_agg(T+1)-agg_ss.K_agg;
   err_tran = err_tran_vec;
   C_init = vpa(0.5*(C_h+C_l));
   iter = iter+1;
   
   if disp_tran==1
       fprintf("iter trans = %d \n",iter)
       fprintf("err  trans = %f \n",err_tran)
       disp(err_tran_vec(end))
       fprintf("C_l = %.18f \n",C_l)
       fprintf("C_h = %.18f \n",C_h)
       toc
       disp('--------------------------------------------')
   end
    
end %END WHILE

% Include dynamic b grid among outputs
agg_tran.b_grid = pol_tran.b_grid;

conv_flag = 0;
if err_tran>tol_tran
    conv_flag = -1;
end

end %END FUNCTION <fun_transition>

