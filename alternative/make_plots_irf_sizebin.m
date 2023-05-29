%% Load results
%% TODO: this script is not up to date - need to add impacted vs unimpacted loops
% Baseline grant
load(fullfile(ResultsDir,'grant_baseline.mat'))
mu_active_baseline = distrib_tran.mu_active;
w_baseline = path.w;
% No grant
load(fullfile(ResultsDir,'nogrant.mat'))
mu_active_nogrant = distrib_tran.mu_active;
w_nogrant = path.w;

% Unpack some parameters
nx = par.nx;
nk = par.nk;
nb = par.nb;
T  = par.T;
ns = par.ns;
A_small = par.A_small;
k_grid  = par.k_grid;
x_grid  = par.x_grid;
eta     = par.eta;
% Steady state wage
wage_ss = sol.wage;
% Define weigths for grant (s=1) and  no grant (s=2) firms
weights= par.weights;

% Define employment cutoffs
bin_lids = [3]';
n_bins = length(bin_lids)+1;
%% Compute Employment share by firm size bins for each t on the transition path
% Distribution: distrib_tran.mu_active and eta
% labor demand by x,kappa,t,grant/nogrant is fun_l(x,wage,kappa,par)
% where x is x_grid(x_c)*A_small(t,i_c)
%       wage is path.w(t)
%       kappa is k_grid(k_c)

emp_sizebins_baseline = zeros(n_bins,T+1);
emp_sizebins_nogrant = zeros(n_bins,T+1);

for i_c = 1:ns  % grant or no grant
    for t   = 1:T+1 % period on the transition path
        wage_baseline = w_baseline(t);
        wage_nogrant = w_nogrant(t);
        for k_c = 1:nk % current capital
            kappa = k_grid(k_c);
            for x_c = 1:nx % current x
                xval = x_grid(x_c)*A_small(t,i_c);
                % Compute firm size given x_c, k_c,t,i_c
                lval_baseline = fun.fun_l(xval,wage_baseline,kappa,par);
                lval_nogrant = fun.fun_l(xval,wage_nogrant,kappa,par);
                % Bin firms according to their size in the steady state!
                lval_ss = fun.fun_l(x_grid(x_c),wage_ss,kappa,par);
                for b_c = 1:nb % current debt
                    % Update emp_sizebins for each policy environment
                    if lval_ss<bin_lids(1)
                        %less than 50 workers in the steady state!
                        emp_sizebins_baseline(1,t) = emp_sizebins_baseline(1,t) + ...
                            lval_baseline*mu_active_baseline(b_c,x_c,k_c,t,i_c)*weights(x_c,k_c,i_c); 
                        emp_sizebins_nogrant(1,t) = emp_sizebins_nogrant(1,t) + ...
                            lval_nogrant*mu_active_nogrant(b_c,x_c,k_c,t,i_c)*weights(x_c,k_c,i_c); 
                    else
                         % 50 or more workers in the steady state!
                        emp_sizebins_baseline(2,t) = emp_sizebins_baseline(2,t) + ...
                            lval_baseline*mu_active_baseline(b_c,x_c,k_c,t,i_c)*weights(x_c,k_c,i_c);
                         emp_sizebins_nogrant(2,t) = emp_sizebins_nogrant(2,t) + ...
                            lval_nogrant*mu_active_nogrant(b_c,x_c,k_c,t,i_c)*weights(x_c,k_c,i_c); 
                    end
                
                end 
            end
        end
    end
end

%% Compute Employment share by firm size bins in the steady state


% Compute employment shre by firm size bin in the steady state
emp_sizebins_ss = zeros(2,1);
for k_c = 1:nk % current capital
    kappa = k_grid(k_c);
    for x_c = 1:nx % current x
        xval = x_grid(x_c);
        for b_c = 1:nb % current debt
            lval = fun.fun_l(xval,wage_ss,kappa,par);
            if lval<bin_lids(1)
                %less than 50 workers
                emp_sizebins_ss(1)   =emp_sizebins_ss(1)+lval*distribS.mu_active(b_c,x_c,k_c);         
            else
                % 50 or more workers
                emp_sizebins_ss(2)   =emp_sizebins_ss(2)+lval*distribS.mu_active(b_c,x_c,k_c);
            end
        end
    end
end

%% Compute impulse response of employment by firm size bins
for bin_c = 1:n_bins
    irf_emp_sizebins_baseline(bin_c,:)= emp_sizebins_baseline(bin_c,:)/emp_sizebins_ss(bin_c);
    irf_emp_sizebins_nogrant(bin_c,:)= emp_sizebins_nogrant(bin_c,:)/emp_sizebins_ss(bin_c);

end
figure
plot(irf_emp_sizebins_baseline(:,1:32)')
figure
plot(irf_emp_sizebins_nogrant(:,1:32)')


