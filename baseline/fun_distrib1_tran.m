function [distrib_tran] = fun_distrib1_tran(par,pol_tran,distrib_ss,sol_ss)
% Purpose: compute the stationary distributions mu_0(k,b,x) before exit
% and the distribution of active firms, mu_active, along the transition
% path. Forward iteration starting from the steady-state mu_0 and mu_active.
%------------------------- input checks ----------------------------------%
if isstruct(par)==0
    error("Input <par> in fun_distrib1_tran must be a structure!")
end
if isstruct(pol_tran)==0
    error("Input <pol_tran> in fun_distrib1_tran must be a structure!")
end
if isstruct(distrib_ss)==0
    error("Input <distrib_ss> in fun_distrib1_tran must be a structure!")
end
%-------------------------------------------------------------------------%

% Unpack policy functions for all t
pol_debt   = pol_tran.pol_debt;   %dim:(nk,nb,nx,T+1,nn)
pol_exit   = pol_tran.pol_exit;   %dim:(nk,nb,nx,T+1,nn)
pol_entry  = pol_tran.pol_entry;  %dim:(nk,nb,nx,T+1,nn)
pol_kp_ind = pol_tran.pol_kp_ind; %dim:(nk,nb,nx,T+1,nn)
b_grid     = pol_tran.b_grid;    %dim:(nk,nb,T+1,nn)
phi_dist   = pol_tran.phi_dist;  %dim:(nk,nb,nx,T+1,nn)

% Unpack the steady-state distribution --> initial condition
% struct contains = {mu,mu_active,entry_vec}
mu_init        = distrib_ss.mu;

validateattributes(mu_init, {'double'}, {'finite', 'nonnan', 'nonempty','>=',0,'real','size', [par.nk,par.nb,par.nx]})

% Unpack PAR:
T        = par.T;
nx       = par.nx;
nb       = par.nb;
nk       = par.nk;
nn       = par.nn;   % nn = ni x ns
pi_x     = par.pi_x; % transition matrix for x
mass_vec     = par.mass_vec; % mass of potential entrants, dim (T+1,1)
psi      = par.psi;
verbose  = par.verbose;

if ~isequal(size(phi_dist),[nk,nb,nx,T+1,nn])
    error('Size of phi_dist not correct')
end

%% Compute mu^0, distrib before entry and exit

mu = zeros(nk,nb,nx,T+1,nn); %dim: (b,x,k,time,impact x grant)

% Initial condition: distribution of incumbents in the first period is the
% steady state distribution
for n_c = 1:nn % impacted vs no impacted x grant vs no grant
    mu(:,:,:,1,n_c) = mu_init;
end

if verbose>=1; disp('Start distribution..'); end

for n_c = 1:nn
    for t=2:T+1
        % mass of potential entrants at t
        mass = mass_vec(t);
        % Policy function at time t to be passed to the subfunction
        pol_debt_t   = pol_debt(:,:,:,t-1,n_c);  %(k,b,x)
        pol_exit_t   = pol_exit(:,:,:,t-1,n_c);  %(k,b,x)
        pol_entry_t  = pol_entry(:,:,:,t-1,n_c); %(k,b,x)
        pol_kp_ind_t = pol_kp_ind(:,:,:,t-1,n_c); %(k,b,x)
        phi_dist_t   = phi_dist(:,:,:,t-1,n_c);  %(k,b,x)
        b_grid_t     = b_grid(:,:,t-1,n_c);      %(k,b)
        
        b_min_all = min(b_grid_t(:,1));   % scalar
        b_max_all = max(b_grid_t(:,nb));  % scalar
        b_grid_all = linspace(b_min_all,b_max_all,nb)';
        
        bopt  = pol_debt_t;                    % dim is (nk,nb,nx)
        b_min = reshape(b_grid_t(pol_kp_ind_t,1),[nk,nb,nx]);
        b_max = reshape(b_grid_t(pol_kp_ind_t,nb),[nk,nb,nx]);
        
        bopt_new = (bopt-b_min)./(b_max-b_min).*(b_max_all-b_min_all)+b_min_all; % dim is (nk,nb,nx)
        
        [loc_vec,omega_vec] = find_loc_vec(b_grid_all,bopt_new(:)); % dim is (nk*nb*nx,1)
        
        left_loc_arr  = reshape(loc_vec,[nk,nb,nx]);
        omega_arr     = reshape(omega_vec,[nk,nb,nx]);
        
        % Iterate forward on distribution equation
        mu(:,:,:,t,n_c) = sub_mu_onestep(mu(:,:,:,t-1,n_c),phi_dist_t,pol_kp_ind_t,...
            pol_exit_t,pol_entry_t,left_loc_arr,omega_arr,pi_x,mass,psi);
        
        if verbose>=2
            fprintf('iter = %d ; n_c = %d \n', t,n_c);
        end
    end %for t
end % for n_c


validateattributes(mu, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0,'size', [nk,nb,nx,T+1,nn]})

%% Compute mu, active firms

% Measure of active firms (see eq. 5)
mu_active = zeros(nk,nb,nx,T+1,nn); %dim:(k,b,x,time,impact x grant)

for n_c = 1:nn
    for t = 1:T+1
        for x_c = 1:nx
            for b_c = 1:nb
                for k_c = 1:nk
                    mu_active(k_c,b_c,x_c,t,n_c)=(1-psi)*(1-pol_exit(k_c,b_c,x_c,t,n_c))*mu(k_c,b_c,x_c,t,n_c) ...
                        + mass*pol_entry(k_c,b_c,x_c,t,n_c)*phi_dist(k_c,b_c,x_c,t,n_c);
                end
            end
        end
    end
end

validateattributes(mu_active, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0,'size', [nk,nb,nx,T+1,nn]})

%% Pack outputs into a struct

distrib_tran = v2struct(mu,mu_active);

end %END FUNCTION <fun_distrib1_tran>
