%% Computation of forced exit in the steady-state and in the transition
clear
clc
close all
% Add path to numerical tools
addpath(genpath(fullfile('tools')));

% %% Exit in steady-state
% disp('Compute fraction of forced exit in the steady-state')
% % Load results
% load(fullfile('mat','ss.mat'),'agg','b_grid','sol','par','prices','distribS','model_mom','data_mom')
% 
% nx = par.nx;
% nb = par.nb;
% nk = par.nk;
% mu = distribS.mu;
% pol_exit_vol    = sol.pol_exit_vol;
% pol_exit_forced = sol.pol_exit_forced;
% 
% % Fraction of total exit due to forced exit
% exit_forced = 0;
% exit_vol = 0;
% for x_c = 1:nx % current x
%     for b_c = 1:nb % current debt
%         for k_c = 1:nk % current capital
%             % forced exit
%             dexit_forced = pol_exit_forced(k_c,b_c,x_c);
%             exit_forced = exit_forced+dexit_forced*mu(k_c,b_c,x_c);
%             % voluntary exit
%             dexit_vol = pol_exit_vol(k_c,b_c,x_c);
%             exit_vol = exit_vol+dexit_vol*mu(k_c,b_c,x_c);
%         end
%     end
% end
% exit_forced = exit_forced/sum(mu(:));
% exit_vol    = exit_vol/sum(mu(:));
% frac_exit_forced = exit_forced/agg.exit_rate;
% frac_exit_vol = exit_vol/agg.exit_rate;
% 
% fprintf('Fraction of forced exit = %f \n',frac_exit_forced)
% fprintf('Fraction of volunt. exit = %f \n',frac_exit_vol)

%% Exit in the transition
clear
clc
disp('Compute fraction of forced exit in the transition')
% Load results
load(fullfile('mat','grant_baseline.mat'),'par','agg_tran','prices','path','pol_tran','distrib_tran','val_tran')
%val_tran is 5-D (k,b,x,1:5,nn)

nx = par.nx;
nb = par.nb;
nk = par.nk;
nn = par.nn;
ni = par.ni;
ns = par.ns;
Xp = par.Xp;
k_grid = par.k_grid;
x_grid = par.x_grid;
A_small = par.A_small;
fixcost = par.fixcost;
wage_ss = prices.wage;
theta = par.theta;
delta = par.delta_k;
b_grid = pol_tran.b_grid; %(k,b,t,n)
mu = distrib_tran.mu; %(k,b,x,t,n)
weights = par.weights; %(k,x,nn)
% val_tran is (k,b,x,1:5,nn)
T_short = size(val_tran,4);
% Compute pol_exit_forced,pol_exit_vol
pol_exit_vol    = zeros(nk,nb,nx,T_short,nn);
pol_exit_forced = zeros(nk,nb,nx,T_short,nn);

% Compute grant
grant_vec  = zeros(nk,nx,ns); %dim: (nk,nx,ns)
for x_c = 1:nx
    for k_c = 1:nk
        kappa = k_grid(k_c);
        % Grant amount is based on labor demand *if there were NO TFP
        % shock (so, doesn't depend on impactedness)
        x_val = x_grid(x_c);
        grant_vec(k_c,x_c,1) = Xp*wage_ss*fun.fun_l(x_val,wage_ss,kappa,par);
    end
end

% Compute policy function for forced vs voluntary exits
for n_c = 1:nn
    [i_c,s_c] = ind2sub([ni,ns],n_c); % i_c = impact indicator; s_c = grant indicator

    for t_c = 1:T_short
        for x_c = 1:nx
            x_val = x_grid(x_c)*A_small(t_c,i_c);
            for b_c = 1:nb
                for k_c = 1:nk
                    kappa = k_grid(k_c);
                    c     = fixcost(k_c);    
                    b_val = b_grid(k_c,b_c,t_c,n_c);
                    % Compute profit
                    profit_val = fun.fun_profit(x_val,kappa,c,path.w(t_c),par);
                    % Add grant to profit
                    if t_c==1
                        profit_val = profit_val+grant_vec(k_c,x_c,s_c);
                    end
                    % Forced liquidation
                    aux1 = profit_val-b_val+theta*(1-delta)*kappa<0;
                    pol_exit_forced(k_c,b_c,x_c,t_c,n_c) = double(aux1);
                    % Voluntary liquidation
                    aux2 = val_tran(k_c,b_c,x_c,t_c,n_c)< theta*(1-delta)*kappa-b_val;
                    pol_exit_vol(k_c,b_c,x_c,t_c,n_c) = double(aux2);
                   
                end
            end
        end
    end
end


% Fraction of total exit due to forced exit
exit_forced = zeros(T_short,1);
exit_vol = zeros(T_short,1);

for t_c = 1:T_short
    exit_forced(t_c) = 0;
    exit_vol(t_c)    = 0;
    for n_c = 1:nn
        for x_c = 1:nx % current x
            for b_c = 1:nb % current debt
                for k_c = 1:nk % current capital
                    % forced exit
                    dexit_forced = pol_exit_forced(k_c,b_c,x_c,t_c,n_c);
                    exit_forced(t_c) = exit_forced(t_c)+dexit_forced*weights(k_c,x_c,n_c)*mu(k_c,b_c,x_c,t_c,n_c);
                    % voluntary exit
                    if dexit_forced==0
                        dexit_vol = pol_exit_vol(k_c,b_c,x_c,t_c,n_c);
                    else
                        dexit_vol=0;
                    end
                    exit_vol(t_c) = exit_vol(t_c)+dexit_vol*weights(k_c,x_c,n_c)*mu(k_c,b_c,x_c,t_c,n_c);
                end
            end
        end
    end
end

% Dimension: (T_short,1)
frac_exit_forced = exit_forced./agg_tran.exit_vec(1:T_short,1);
frac_exit_vol = exit_vol./agg_tran.exit_vec(1:T_short,1);


