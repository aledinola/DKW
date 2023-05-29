function [mom] = fun_targets(sol,mustruct,par,prices,agg,b_grid)
% Purpose: 
% fun_targets computes calibration targets.
% #VC# V53
% INPUTS
% sol      : struct with solution policies from vfi
% mustruct : struct with distributions
% par      : struct with model parameters
% prices   : struct with prices
% agg      : struct with aggregate moments
% b_grid   : Grid for debt, dim (nk,nb)
% OUTPUTS
% mom      : Struct with model moments

%------------------------- input checks ----------------------------------%
if ~isstruct(sol)
    error("Input sol must be a structure")
end
if ~isstruct(mustruct)
    error("Input mustruct in fmust be a structure")
end
if ~isstruct(par)
    error("Input par in must be a structure")
end
if ~isstruct(prices)
    error("Input prices in must be a structure")
end
if ~isstruct(agg)
    error("Input agg in must be a structure")
end
%-------------------------------------------------------------------------%

mom = struct();

% Unpack relevant variables:
pol_debt   = sol.pol_debt; %dim: (nk,nb,nx)
pol_exit   = sol.pol_exit;
pol_kp_ind = sol.pol_kp_ind;
%pol_exit_vol    = sol.pol_exit_vol;
pol_exit_forced = sol.pol_exit_forced;
pol_exit_vol    = sol.pol_exit_vol;

pol_entry = sol.pol_entry;
phi_dist  = sol.phi_dist;
mu        = mustruct.mu;
mu_active = mustruct.mu_active;
wage      = prices.wage;

% Unpack some parameters:
pi_x   = par.pi_x;
nx     = par.nx;
nb     = par.nb;
nk     = par.nk;
x_grid = par.x_grid;
k_grid = par.k_grid;
mass   = par.mass;
psi    = par.psi;
fixcost = par.fixcost;
N_sim = par.N_sim;
T_sim = par.T_sim;

if ~isequal(size(b_grid),[nk,nb])
    error('Size of b_grid is wrong')
end

%% Define total exit rate and usefule permutations

exit_all = psi+(1-psi)*pol_exit;

% NOTE: suffix "p" stands for "permuted"

% exit_all(k,b,x) ==> exit_allp(b,x,k)
exit_allp = permute(exit_all,[2 3 1]);
% pol_debt(k,b,x) ==> pol_debt(b,x,k)
%pol_debtp = permute(pol_debt,[2 3 1]);
% b_grid(k,b) ==> b_grid(b,k)
b_gridp = permute(b_grid,[2 1]);
% labor_vec(k,b) ==> labor_vec(b,k)
% labor_vecp = permute(labor_vec,[2 1]);

%% Output or revenue share small firms 
% TODO: do we include the fixcost in revenue small?
revshare_small = agg.output_small/(agg.output_small+agg.Y_corp);

%% Employment share small firms
empshare_small = (agg.L_agg-agg.L_corp)/(agg.L_agg);

%% Share of firms with positive debt

b_mat = repmat(b_grid,[1,1,par.nx]);
hasNetDebt = sum(mu_active(b_mat>0),"all") / sum(mu_active,'all');

%% Precompute labor demand on the x grid

labor_vec   = zeros(nk,nx);
revenue_vec = zeros(nk,nx);
for x_c = 1:nx
    x_val = x_grid(x_c);
    for k_c = 1:nk
        kappa = k_grid(k_c);
        c     = fixcost(k_c); 
        labor_vec(k_c,x_c)  = fun.fun_l(x_val,wage,kappa,par);
        revenue_vec(k_c,x_c)= fun.prod_small(x_val,kappa,labor_vec(k_c,x_c),c,par);
    end
end

% Employment at active firms
empl_active = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            empl_active = empl_active+labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
        end
    end
end

% labor_vec(k,b) ==> labor_vec(b,k)
labor_vecp = permute(labor_vec,[2 1]);

%% Average firm size (small firms)

avefirmsize = empl_active/sum(mu_active(:));

%% Average firm size by firm age

% Age 0
empl_age0 = 0;

for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            empl_age0=empl_age0+mass*labor_vec(k_c,x_c)*pol_entry(k_c,b_c,x_c)*phi_dist(k_c,b_c,x_c);
        end
    end
end

% We already computed this as Mentr in fun_obj.m
mass_entrants = mass*sum(pol_entry.*phi_dist,'all');

avefirmsize_age0 = empl_age0/mass_entrants;

%% Job creation and job destruction of continuing firms, autocorrelation of empl 
% b_gridp has dim (nb,nk)
% labor_vecp has dim (nx,nk)

% Precompute interpolated exit policy function dexit(x',k,b,x)
% stay_arr1 = zeros(nx,nk,nb,nx);
% for x_c = 1:nx % current x
%     for b_c = 1:nb % current debt
%         for k_c = 1:nk % current capital
%             bnext = pol_debt(k_c,b_c,x_c);
%             knext_ind = pol_kp_ind(k_c,b_c,x_c);
%             exit_all_bx = exit_allp(:,:,knext_ind); % dim: (nb,nx)
%             dexit_inter = myinterp1q(b_gridp(:,knext_ind),exit_all_bx,bnext); % dim is (1,nx')
%             dexit = min(max(dexit_inter,0),1); % dim is (1,nx')
%             stay_arr1(:,k_c,b_c,x_c) = 1-dexit;
%         end %k_c
%     end %b_c
% end %x_c

% Faster code. Assumption: grid for b must be equally spaced
K = pol_kp_ind;
bminK = reshape(b_gridp(1,K),size(K));
bmaxK = reshape(b_gridp(nb,K),size(K));
Y = 1 + (nb-1) * (pol_debt - bminK) ./ (bmaxK-bminK);
Yt = max(min(Y,nb-1),1); % no need if there is no overflowed in the data
I = floor(Yt); % nk x nb x nx
W = Y-I;
[I,J]=ndgrid(I,1:nx); % (nk x nb x nx) x nx
K = reshape(repmat(K,[1 1 1 nx]),size(I));
rhsilin = sub2ind(size(exit_allp),I,J,K); % (nk x nb x nx) x nx;
rhsilin = reshape(rhsilin, [nk,nb,nx,nx]);
dexit_inter = (1-W).*exit_allp(rhsilin) + W.*exit_allp(rhsilin+1);
dexit_inter = permute(dexit_inter, [4 1 2 3]); % [nx,nk,nb,nx]
dexit = min(max(dexit_inter,0),1);
stay_arr = 1-dexit;

%err = max(abs(stay_arr-stay_arr1),[],'all')

jc_num       = 0; % Numerator in the definition of JC rate
jd_num       = 0; % Numerator in the definition of JD rate
mean_l       = 0; % Average of employment in surviving firms
contin_firms = 0; % Mass of continuing firms

for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            knext_ind = pol_kp_ind(k_c,b_c,x_c);
            stay = stay_arr(:,k_c,b_c,x_c)';% dim is (1,nx')
            
            % Job creation
            new_jobs = max(0,labor_vecp(:,knext_ind)-labor_vec(k_c,x_c));% (nx',1)
            temp = pi_x(x_c,:).*new_jobs'.*stay; % (1,nx')
            jc_num = jc_num+sum(temp)*mu_active(k_c,b_c,x_c);
            
            % Job destruction
            firing = max(0,labor_vec(k_c,x_c)-labor_vecp(:,knext_ind));% (nx',1)
            temp = pi_x(x_c,:).*firing'.*stay; % dim is (1,nx')
            jd_num = jd_num+sum(temp)*mu_active(k_c,b_c,x_c);
            
            % Average of employment in surviving firms
            temp = labor_vecp(:,knext_ind)'.*stay.*pi_x(x_c,:); % dim is (1,nx')
            mean_l = mean_l + sum(temp)*mu_active(k_c,b_c,x_c);
            temp2 = stay.*pi_x(x_c,:); % dim is (1,nx')
            contin_firms = contin_firms+sum(temp2)*mu_active(k_c,b_c,x_c);  
        end %k_c
    end %b_c
end %x_c

jc_rate = jc_num/empl_active;
jd_rate = jd_num/empl_active;
% NOTE: mean_l should be normalized by the mass of continuing firms
mean_l  = mean_l/contin_firms;

var_l        = 0; % Variance of employment in surviving firms
cov_l        = 0; % Covariance l,l' in surviving firms

for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            knext_ind = pol_kp_ind(k_c,b_c,x_c);
            stay = stay_arr(:,k_c,b_c,x_c)';% dim is (1,nx')
            
            % Variance of employment in surviving firms, covariance l,l' in surviving firms
            var_l = var_l + sum((labor_vecp(:,knext_ind)'-mean_l).^2.*stay.*pi_x(x_c,:))*mu_active(k_c,b_c,x_c); %scalar
            temp = (labor_vecp(:,knext_ind)'-mean_l).*stay.*pi_x(x_c,:); % dim is (1,nx')
            cov_l = cov_l + (labor_vec(k_c,x_c)-mean_l)*sum(temp)*mu_active(k_c,b_c,x_c); %scalar
            
        end %k_c
    end %b_c
end %x_c

% Correlation l,l' in surviving firms
corr_l = cov_l/var_l;

%% Compute firm size (employment) distribution

% We want 
% share of firms 
% share of employment 
% for these size classes:
% [0,5), [5, 10), [10, 20), >=20
% TODO: change classes to 
% [0,5), [5, 10), [10, 20), [20, 100), >=100


bin_lids = [10, 20, 100]';

firms    = zeros(length(bin_lids)+1,1);
workers  = zeros(length(bin_lids)+1,1);
revenues = zeros(length(bin_lids)+1,1);

% Vectorized code **********************************
%[share_firms_bis,share_empl_bis,share_rev_bis] = fun_bins(mu_active,labor_vec,revenue_vec,bin_lids);
% * End vectorized code ***************************

for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            lval = labor_vec(k_c,x_c);
            
            if lval<bin_lids(1)
                %less than 10 workers
                firms(1)   =firms(1)+mu_active(k_c,b_c,x_c);
                workers(1) =workers(1)+labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
                revenues(1)=revenues(1)+revenue_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
                
            elseif lval>=bin_lids(1) && lval<bin_lids(2)
                % [10,20)
                firms(2)   = firms(2)+mu_active(k_c,b_c,x_c);
                workers(2) = workers(2)+labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
                revenues(2)= revenues(2)+revenue_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
                
            elseif lval>=bin_lids(2) && lval<bin_lids(3)
                % [20,100)
                firms(3)    = firms(3)+mu_active(k_c,b_c,x_c);
                workers(3)  = workers(3)+labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
                revenues(3) = revenues(3)+revenue_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
                                
            elseif lval>=bin_lids(3)
                % >= 100
                firms(4)    = firms(4)+mu_active(k_c,b_c,x_c);
                workers(4)  = workers(4)+labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
                revenues(4) = revenues(4)+revenue_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
            else
                error("smth wrong")
            end
        end
    end
end

share_firms = firms/sum(firms);
share_empl  = workers/sum(workers);
share_rev   = revenues/sum(revenues);

%% Exit and entry rates

% Fraction of total exit due to forced exit
exit_forced = 0;
exit_vol = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            % forced exit
            dexit_forced = pol_exit_forced(k_c,b_c,x_c);
            exit_forced = exit_forced+dexit_forced*mu(k_c,b_c,x_c);
            % voluntary exit
            dexit_vol = pol_exit_vol(k_c,b_c,x_c);
            exit_vol = exit_vol+dexit_vol*mu(k_c,b_c,x_c);
        end
    end
end
exit_forced = exit_forced/sum(mu(:));
exit_vol    = exit_vol/sum(mu(:));
frac_exit_forced = exit_forced/agg.exit_rate;
frac_exit_vol = exit_vol/agg.exit_rate;

entry_rate = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            entry = mass*pol_entry(k_c,b_c,x_c);
            entry_rate = entry_rate+entry*phi_dist(k_c,b_c,x_c);
        end
    end
end
entry_rate = entry_rate/sum(mu(:));

%% Exit rate by firm size
bin_lids = [10, 20, 100]';
%bin_lids = [0.1, 1, 10]';
min_size = par.emp_min; % we exclude microbusinesses
%min_size=0;
exit_rate_size = zeros(length(bin_lids)+1,1);
mass_size      = zeros(length(bin_lids)+1,1);

for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            lval  = labor_vec(k_c,x_c);
            dexit = exit_all(k_c,b_c,x_c);
            
            if lval>min_size && lval<bin_lids(1)
                %less than 10 workers (but more than min_size)
                
                exit_rate_size(1) =  exit_rate_size(1)+dexit*mu(k_c,b_c,x_c);
                mass_size(1) = mass_size(1) + mu(k_c,b_c,x_c);
            elseif lval>=bin_lids(1) && lval<bin_lids(2)
                % [10,20)
                exit_rate_size(2) =  exit_rate_size(2)+dexit*mu(k_c,b_c,x_c);
                mass_size(2) = mass_size(2) + mu(k_c,b_c,x_c);
            elseif lval>=bin_lids(2) && lval<bin_lids(3)
                % [20,100)
                exit_rate_size(3) =  exit_rate_size(3)+dexit*mu(k_c,b_c,x_c);
                mass_size(3) = mass_size(3) + mu(k_c,b_c,x_c);
            elseif lval>=bin_lids(3)
                % 100+
                exit_rate_size(4) =  exit_rate_size(4)+dexit*mu(k_c,b_c,x_c);
                mass_size(4) = mass_size(4) + mu(k_c,b_c,x_c);
            end
            
        end
    end
end

exit_rate_size = exit_rate_size./(1e-20+mass_size);


%% Leverage of small firms (defined as debt/assets)

% If firm has debt, i.e. b>0, leverage is b/owned capital
% If firm has financial assets, i.e. b<0, leverage is zero

leverage = 0;
assets   = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            kappa = k_grid(k_c);
            leverage = leverage+max(b_grid(k_c,b_c),0)*mu_active(k_c,b_c,x_c);
            assets   = assets + kappa*mu_active(k_c,b_c,x_c);
        end
    end
end
leverage = leverage/assets;

% Leverage for entrants
leverage_entrants = 0;
assets_entrants   = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            kappa = k_grid(k_c);
            leverage_entrants = leverage_entrants+max(b_grid(k_c,b_c),0)*mass*pol_entry(k_c,b_c,x_c)*phi_dist(k_c,b_c,x_c);
            assets_entrants   = assets_entrants+kappa*mass*pol_entry(k_c,b_c,x_c)*phi_dist(k_c,b_c,x_c);
        end
    end
end
leverage_entrants = leverage_entrants/assets_entrants;

%% Labor share small firms
% TODO: should we include fixed cost in output (denominator of labor share)
payroll = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            payroll = payroll+wage*labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
        end
    end
end

laborshare_small = payroll/agg.output_small;

%%  Debt to payroll ratio, given positive debt

tot_debt_pos = 0;
payroll_debt = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            b_val = b_grid(k_c,b_c);
            if b_val>0
                tot_debt_pos = tot_debt_pos+b_val*mu_active(k_c,b_c,x_c);
                payroll_debt = payroll_debt +wage*labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
            end
        end
    end
end
if payroll_debt<=0
    debt_payroll_cond=0;
else
    debt_payroll_cond = tot_debt_pos/payroll_debt;
end

%%  Cash to payroll ratio, conditional on having cash
tot_debt_neg = 0;
payroll_cash = 0;
for x_c = 1:nx % current x
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            b_val = b_grid(k_c,b_c);
            if b_val<=0
                tot_debt_neg = tot_debt_neg+abs(b_val)*mu_active(k_c,b_c,x_c);
                payroll_cash = payroll_cash +wage*labor_vec(k_c,x_c)*mu_active(k_c,b_c,x_c);
            end
        end
    end
end
if (payroll_cash<=0) 
    cash_payroll_cond = 0;
else
    cash_payroll_cond = tot_debt_neg/payroll_cash;
end

%% Capital-to-revenue ratio for small firms

capital = 0;
for k_c = 1:nk % current capital
   kappa = k_grid(k_c);
   capital = capital+kappa*sum(mu_active(k_c,:,:),'all');
end
tot_fixcost = 0;
% Compute total fixed cost
for x_c = 1:nx % current capital
    tot_fixcost = tot_fixcost+fixcost(k_c)*sum(mu_active(:,:,x_c),'all');
end
% Small firm revnue = production output
temp = tot_fixcost+agg.output_small;
caprev = capital/temp;

%% Capital-to-payroll ratio for small firms
k_payroll_ratio = capital/payroll;

%% Fixed cost to payroll, ratio
% Total fixed cost: tot_fixcost
fixedcost_to_payroll = tot_fixcost/payroll;

%% Fixed cost to revenues, ratio
% Note that revenue = production outputs (need to add back fixed costs to
% output_small since fixed costs are subtracted).
fixedcost_to_rev = tot_fixcost/(tot_fixcost+agg.output_small);

%% Investment rate of unconstrained firms
% Indicator for unconstrained firms: 
% Unconstrained firms are those such that b <= B_hat(k,x)
B_hat_mat = repmat(sol.B_hat,[1,1,nb]);
B_hat_mat = permute(B_hat_mat,[1,3,2]);
is_uc = (b_mat<B_hat_mat); % (nk,nb,nx)

sim = fun_simulate(sol,par,mustruct,is_uc,b_grid,N_sim,T_sim);

% k_sim is an index, k_sim_val=k_grid(k_sim)
k_sim = sim.k_sim_val;  
%b_sim = sim.b_sim;
%x_sim = sim.x_sim;
surv_sim = sim.surv_sim;

% Keep only firms that survives
k_sim = k_sim(surv_sim==1,4:4:T_sim);
%b_sim = b_sim(surv_sim==1,4:4:T_sim);
%x_sim = x_sim(surv_sim==1,4:4:T_sim);

% Number of surviving firms
N_surv = size(k_sim,1);
% Number of simulation years
Y_sim = size(k_sim,2);

% Investment rate (k'-(1-delta)k)/k
invrate = (k_sim(:,2:Y_sim)-(1-par.delta_k)*k_sim(:,1:Y_sim-1))./k_sim(:,1:Y_sim-1);

% mean and standard deviation of invrate
mean_invrate  = sum(invrate,'all')/numel(invrate);
stdev_invrate = sum(std(invrate,0,2))/N_surv; 

% serial correlation of investment rate
scor_invrate_vec = zeros(N_surv,1);
for i_c = 1:N_surv
    % corrcoef gives NaN sometimes
    var_invrate = cov(invrate(i_c,2:end),invrate(i_c,1:end-1));
    scor_invrate_vec(i_c) = var_invrate(2,1)/ (var_invrate(1,1)+1e-20);
    %R_aux = corrcoef(invrate(i_c,2:end),invrate(i_c,1:end-1));
    %scor_invrate_vec(i_c) = R_aux(1,2);
end
scor_invrate = sum(scor_invrate_vec)/numel(scor_invrate_vec);

% Frequency of lumpy investment or positive inv spike (investment rate>=20%)
lumpy_inv = invrate>0.2;
freq_lumpinv   = sum(lumpy_inv,'all')/numel(lumpy_inv);

%% Pack moments into struc
mom.avefirmsize          = avefirmsize;
mom.empshare_small       = empshare_small;
mom.revshare_small       = revshare_small;
mom.payroll_VA_ratio     = laborshare_small;
%mom.exitrate             = agg.exit_rate;
mom.exitrate             = agg.exit_rate_emp; % exit rate excluding micro-businesses
mom.frac_exit_forced     = frac_exit_forced;
mom.frac_exit_vol     = frac_exit_vol;

mom.exitrate_0_9         = exit_rate_size(1);
mom.exitrate_10_19       = exit_rate_size(2);
mom.exitrate_20_99       = exit_rate_size(3);
mom.exitrate_100_499     = exit_rate_size(4);

mom.entryrate            = entry_rate;
mom.avefirmsize_age0     = avefirmsize_age0;

mom.jcr                  = jc_rate;
mom.jdr                  = jd_rate;
%mom.jc_rate_entry        = jc_rate_entry;
%mom.jd_rate_exit         = jd_rate_exit;
mom.fixedcost_to_rev     = fixedcost_to_rev;
mom.fixedcost_to_payroll = fixedcost_to_payroll;
mom.autocorr_emp         = corr_l;
mom.debt_asset_all       = leverage; 
mom.debt_asset_entrants  = leverage_entrants; 
mom.debt_asset_all_entrants  = leverage/leverage_entrants; 

mom.k_VA_ratio		     = caprev;
mom.k_payroll_ratio      = k_payroll_ratio;
mom.ave_work             = agg.L_agg;

mom.hasNetDebt = hasNetDebt;
mom.debt_payroll_cond = debt_payroll_cond;
mom.cash_payroll_cond = cash_payroll_cond;


mom.firmshare_0_9    = share_firms(1);
mom.firmshare_10_19   = share_firms(2);
mom.firmshare_20_99   = share_firms(3);
mom.firmshare_100_499 = share_firms(4);     

mom.empshare_0_9     = share_empl(1);
mom.empshare_10_19   = share_empl(2);
mom.empshare_20_99   = share_empl(3);
mom.empshare_100_499 = share_empl(4);

mom.revshare_0_9     = share_rev(1);
mom.revshare_10_19   = share_rev(2);
mom.revshare_20_99   = share_rev(3);
mom.revshare_100_499 = share_rev(4);

% These are 4*1
mom.share_firms = share_firms;
mom.share_empl  = share_empl;
mom.share_rev   = share_rev;

% Investment rate moments
mom.stdev_invrate = stdev_invrate;
mom.scor_invrate = scor_invrate;
mom.freq_lumpinv = freq_lumpinv;
mom.mean_invrate = mean_invrate;

end %end function fun_targets_vec1

