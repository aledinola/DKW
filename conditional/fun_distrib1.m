function [mu,mu_active,entry_vec,flag_mu,dist,iter] = fun_distrib1(par,sol,b_grid,phi_dist)

% Purpose: compute the stationary distribution mu_0(k,b,x), called mu in the code,
% and the distribution of active firms, mu_active(k,b,x).
% Note: b_grid is (nk,nb) array
%       phi_dist is (nk,nb,nx) array
% This function calls:
%   - mu_onestep_mex
%------------------------- input checks ----------------------------------%
if isstruct(par)==0
    error("Input <par> must be a structure!")
end
if isstruct(sol)==0 
    error("Input <sol> must be a structure!")
end
% validateattributes(phi_dist, {'double'}, {'finite', 'nonnan', 'nonempty',...
%     'real', '>=', 0,'size', [par.nb,par.nx]})
%-------------------------------------------------------------------------%
 
% Unpack SOL:
pol_kp_ind = sol.pol_kp_ind;% dim: (nk,nb,nx)
pol_debt   = sol.pol_debt;  % dim: (nk,nb,nx)
pol_exit   = sol.pol_exit;  % dim: (nk,nb,nx)
pol_entry  = sol.pol_entry; % dim: (nk,nb,nx)

validateattributes(pol_kp_ind, {'double'}, {'finite', 'nonnan', 'nonempty','real','>=',0,'<=',par.nk,'size', [par.nk,par.nb,par.nx]})
validateattributes(pol_debt, {'double'}, {'finite', 'nonnan', 'nonempty','real','size', [par.nk,par.nb,par.nx]})
validateattributes(pol_exit, {'double'}, {'finite', 'nonnan', 'nonempty','real','>=',0,'<=',1,'size', [par.nk,par.nb,par.nx]})
validateattributes(pol_entry, {'double'}, {'finite', 'nonnan', 'nonempty','real','>=',0,'<=',1,'size', [par.nk,par.nb,par.nx]})
validateattributes(b_grid, {'double'}, {'finite', 'nonnan', 'nonempty','real','size', [par.nk,par.nb]})

% Unpack PAR:
nx       = par.nx;
nb       = par.nb;
nk       = par.nk;
pi_x     = par.pi_x;
mass     = par.mass; %mass of potential entrants
psi      = par.psi;  %exogenous exit rate
tol_dist = par.tol_dist;
maxit    = par.maxiter_dist;
verbose  = par.verbose;
disp_mu  = par.disp_mu;

% Initial condition
mu = zeros(nk,nb,nx);
for k_c = 1:nk
    for b_c = 1:nb
        mu(k_c,b_c,:) = par.x_prob;
    end
end
%mu = repmat(par.x_prob',[nk,nb,1]);
mu = 0.02*mu/sum(mu,'all');
%mu = 0.02*phi_dist/sum(phi_dist(:));
if ~isequal(size(mu),[nk,nb,nx])
    error('mu has wrong size')
end

nn = nk*nb*nx;

dist = 10;
iter = 0;

flag_mu = 0;

if verbose>=1
    fprintf('--------------------------------------------\n')
    fprintf('DISTRIBUTION \n')
    fprintf('--------------------------------------------\n')
end

%% Compute entry vector

entry_vec = zeros(nk,nb,nx);
for x_c = 1:nx % current productivity
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            knext_ind = pol_kp_ind(k_c,b_c,x_c);
            bnext     = pol_debt(k_c,b_c,x_c);
            [left_loc,omega] = find_loc(b_grid(knext_ind,:)',bnext);
            for xp_c = 1:nx
                entry_vec(knext_ind,left_loc,xp_c) =  entry_vec(knext_ind,left_loc,xp_c)+omega*mass*pi_x(x_c,xp_c)*pol_entry(k_c,b_c,x_c)*phi_dist(k_c,b_c,x_c);
                entry_vec(knext_ind,left_loc+1,xp_c)= entry_vec(knext_ind,left_loc+1,xp_c)+(1-omega)*mass*pi_x(x_c,xp_c)*pol_entry(k_c,b_c,x_c)*phi_dist(k_c,b_c,x_c);
            end
        end
    end
end

%% Compute mu^0(b,x)
% Note: here we could use the entry vector already computed above

if verbose>=1
    fprintf(' \n')
    disp('Start distribution..'); 
end

left_loc_arr = ones(nk,nb,nx);
omega_arr    = zeros(nk,nb,nx);

for x_c = 1:nx % current productivity
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            knext_ind = pol_kp_ind(k_c,b_c,x_c);
            bopt      = pol_debt(k_c,b_c,x_c);
            [left_loc,omega] = find_loc(b_grid(knext_ind,:)',bopt);
            omega_arr(k_c,b_c,x_c)    = omega;
            left_loc_arr(k_c,b_c,x_c) = left_loc;
        end
    end
end

%mu1 = zeros(nk,nb,nx);
while dist>tol_dist && iter<=maxit
    iter = iter+1;
    
    % Iterate on distribution equation
    mu1 = sub_mu_onestep(mu,phi_dist,pol_kp_ind,pol_exit,pol_entry,...
        left_loc_arr,omega_arr,pi_x,mass,psi);
    %mu1 = mu_onestep_mex(mu,phi_dist,int32(pol_kp_ind),pol_debt,pol_exit,...
    %    pol_entry,int32(left_loc_arr),omega_arr,b_grid,pi_x,mass,psi);
    
    
    %Compute error
    dist1 = max(abs(mu(:)-mu1(:)));
    dist = max(abs(mu(:)-mu1(:)))/sum(mu,"all")*nn;
    % Update mu
    mu = mu1;
    
    if disp_mu==1
        fprintf('sum(mu1) = %f\n', sum(mu1(:)));
        fprintf('iter = %d, dist_abs = %10.20f, dist_rel = %10.20f \n', iter,dist1,dist);
    end
    
end %end while

validateattributes(mu, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0,'size', [nk,nb,nx]})

%% Compute mu(b,x), active firms

% Measure of active firms (see eq. 5)
mu_active = zeros(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            mu_active(k_c,b_c,x_c)=(1-psi)*(1-pol_exit(k_c,b_c,x_c))*mu(k_c,b_c,x_c) ...
                + mass*pol_entry(k_c,b_c,x_c)*phi_dist(k_c,b_c,x_c);
        end
    end
end

validateattributes(mu_active, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0,'size', [nk,nb,nx]})
validateattributes(entry_vec, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0,'size', [nk,nb,nx]})

if dist>tol_dist
    flag_mu = -1;
else
    if verbose>=1; fprintf('MU (distr.) converged after %d iterations! \n', iter); end
end

end %END FUNCTION <fun_distrib1>
