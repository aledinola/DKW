function [val_c_new,pol_kp_ind] = sub_vfi_onestep(val_c,val0_u,kp_bar,B_hat,profit_mat,...
    k_grid,b_grid,pi_x,theta,q,delta,psi,do_howard,n_howard)

% DESCRIPTION: One step of the Bellman operator T: val_new = T(val)
% This function is called either in the context of inifinite horizon or
% in the backward iteration for the transition.
%
% INTERPOLATION:
% In an earlier version we used bilinear interpolation of v^0 over k' and
% b'. Now we force k' on k_grid and interpolate only over b'.
%
% INPUTS:
% val_c(k,b,x):  value of a constrained firm after exit
% val0_u(k,b,x): value of an unconstrained firm before exit
% kp_bar(k,b,x): upper bound for next-period k', see eq. (25)
% profit_mat(k,x): static profit, dim: (nk,nx)
% B_hat(k,x):      dim(nk,nx)
% k_grid           dim(nk,1),   fixed grid for capital 
% b_grid           dim(nk,nb),  flexible grid for debt   
% pi_z             dim(nx,nx),  transition matrix for x
% theta,q,delta,psi    real(8) parameters
% do_howard,n_howard integer parameters
%
% OUTPUTS:
% val_c_new: dim (nk,nb,nx), real(8)
% pol_kp_ind:  dim (nk,nb,nx), integer
%-------------------------------------------------------------------------%

[nk,nb,nx] = size(val_c);

if nk~=length(k_grid)
    error('nk is not correct')
end
if nb~=size(b_grid,2)
    error('nb is not correct')
end

% Initialize outputs
val_c_new  = zeros(nk,nb,nx);
pol_kp_ind = ones(nk,nb,nx);

% STEP 2 - Impose liquidation, eq. (24)

% Value of constrained firm before exit
val0_c = zeros(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            k_val = k_grid(k_c);
            b_val = b_grid(k_c,b_c);
            profit_val = profit_mat(k_c,x_c);
            liq1 = profit_val-b_val+theta*(1-delta)*k_val<0;
            liq2 = val_c(k_c,b_c,x_c)<theta*(1-delta)*k_val-b_val;
            if liq1 || liq2
                val0_c(k_c,b_c,x_c) = theta*(1-delta)*k_val-b_val;
            else
                val0_c(k_c,b_c,x_c) = val_c(k_c,b_c,x_c);
            end
        end
    end
end

% STEP 3 - Do equation (23)
val0 = zeros(nk,nb,nx);
is_c = ones(nk,nb,nx); % indicator for constrained firms
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            b_val = b_grid(k_c,b_c);
            if b_val<=B_hat(k_c,x_c)
                val0(k_c,b_c,x_c) = val0_u(k_c,b_c,x_c);
                is_c(k_c,b_c,x_c) = 0;
            else
                val0(k_c,b_c,x_c) = val0_c(k_c,b_c,x_c);
            end
        end
    end
end

% STEP 4 - Solve eq. (25)

% Maximize over k' on the grid

for x_c = 1:nx
    %disp(x_c)
    % Compute expected value in eq. (25)
    EVx = zeros(nk,nb); %(k',b')
    for xp_c = 1:nx
        EVx = EVx+val0(:,:,xp_c)*pi_x(x_c,xp_c);
    end
    for b_c = 1:nb
        for k_c = 1:nk
            k_val = k_grid(k_c);
            b_val = b_grid(k_c,b_c);
            profit_val = profit_mat(k_c,x_c);
            kp_ub = kp_bar(k_c,b_c,x_c);
            kp_ub_ind = find(k_grid<=kp_ub, 1, 'last' );

            if is_c(k_c,b_c,x_c)==1 && profit_val-b_val+theta*(1-delta)*k_val>=0

                % Create the RHS of Bellman eq. 25 for each k' \in k_grid
                rhs_vec = repmat(-100000,nk,1);
                %b_prime_vec = zeros(nk,1);
                for kp_c = 1:kp_ub_ind
                    kp_val    = k_grid(kp_c);
                    b_grid_kp = b_grid(kp_c,:)';
                    bprime = max(b_grid_kp(1),(1/q)*(b_val-profit_val+fun.adjcost_scal(kp_val,k_val,theta,delta)));
                    %bprime_vec(kp_c) = bprime';
                    %bprime = (1/q)*(b_val-profit_val+fun.adjcost_scal(kp_val,k_val,theta,delta));
                    v0_int = myinterp1(b_grid_kp,EVx(kp_c,:)',bprime,1);
                    rhs_vec(kp_c) = q*(psi*(theta*(1-delta)*kp_val-bprime)+(1-psi)*v0_int);
                end
                [fval,max_ind]          = max(rhs_vec);
                pol_kp_ind(k_c,b_c,x_c) = max_ind;
                val_c_new(k_c,b_c,x_c)  = fval;
            else
                % firm is unconstrained or forced liquidation: val_c is irrelevant
                val_c_new(k_c,b_c,x_c) = theta*(1-delta)*k_val-b_val;
            end
        end %k
    end %b
end %x

% Now do Howard

% Howard acceleration
if do_howard==1
    val_c_new = fun_howard(val_c_new,pol_kp_ind,val0_u,kp_bar,B_hat,...
        profit_mat,k_grid,b_grid,pi_x,theta,q,delta,psi,n_howard);
end

end %END FUNCTION "sub_vfi_onestep"
%=========================================================================%

function [val_c_new] = fun_howard(val_c,pol_kp_ind,val0_u,kp_bar,B_hat,...
    profit_mat,k_grid,b_grid,pi_x,theta,q,delta,psi,n_howard)

% DESCRIPTION
% This function does one step of the Howard policy improvement algorithm.
% It is called n_howard times by sub_vfi_onestep.
% NOTES
% The input argument "kp_bar" is not used but we keep it for debugging.
%-------------------------------------------------------------------------%

[nk,nb,nx] = size(val_c);

if nk~=length(k_grid)
    error('nk is not correct')
end
if nb~=size(b_grid,2)
    error('nb is not correct')
end

% Initialize output
val_c_new = zeros(nk,nb,nx);

% Precomputations to speed up howard
%allocate(left_arr(nk,nb,nx),omega_arr(nk,nb,nx),bprime_arr(nk,nb,nx))
left_arr   = ones(nk,nb,nx);
omega_arr  = zeros(nk,nb,nx);
bprime_arr = zeros(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            k_val      = k_grid(k_c);
            b_val      = b_grid(k_c,b_c);
            profit_val = profit_mat(k_c,x_c);
            kp_c       = pol_kp_ind(k_c,b_c,x_c);
            kp_val     = k_grid(kp_c);
            b_grid_kp  = b_grid(kp_c,:)';
            bprime = max(b_grid_kp(1),(1/q)*(b_val-profit_val+fun.adjcost_scal(kp_val,k_val,theta,delta)));
            %! It seems that we don't need to store bprime in bprime_arr
            bprime_arr(k_c,b_c,x_c) = bprime;
            [jstar,omega] = myfind_loc( b_grid_kp,bprime);
            left_arr(k_c,b_c,x_c)  = jstar;
            omega_arr(k_c,b_c,x_c) = omega;
        end
    end
end

for h_c = 1:n_howard

% STEP 2 - Impose liquidation, eq. (24)

% Value of constrained firm before exit
val0_c = zeros(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            k_val = k_grid(k_c);
            b_val = b_grid(k_c,b_c);
            profit_val = profit_mat(k_c,x_c);
            liq1 = profit_val-b_val+theta*(1-delta)*k_val<0;
            liq2 = val_c(k_c,b_c,x_c)<theta*(1-delta)*k_val-b_val;
            if liq1 || liq2
                val0_c(k_c,b_c,x_c) = theta*(1-delta)*k_val-b_val;
            else
                val0_c(k_c,b_c,x_c) = val_c(k_c,b_c,x_c);
            end
        end
    end
end

% STEP 3 - Do equation (23)
val0 = zeros(nk,nb,nx);
is_c = ones(nk,nb,nx); % indicator for constrained firms

for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            b_val = b_grid(k_c,b_c);
            if b_val<=B_hat(k_c,x_c)
                val0(k_c,b_c,x_c) = val0_u(k_c,b_c,x_c);
                is_c(k_c,b_c,x_c) = 0;
            else
                val0(k_c,b_c,x_c) = val0_c(k_c,b_c,x_c);
            end
        end
    end
end

% STEP 4 - Solve eq. (25)

for x_c = 1:nx
    %disp(x_c)
    % Compute expected value in eq. (25)
    EVx = zeros(nk,nb); %(k',b')
    for xp_c = 1:nx
        EVx = EVx+pi_x(x_c,xp_c)*val0(:,:,xp_c);
    end
    for b_c = 1:nb
        for k_c = 1:nk
            k_val = k_grid(k_c);
            b_val = b_grid(k_c,b_c);
            profit_val = profit_mat(k_c,x_c);
            %kp_ub = kp_bar(k_c,b_c,x_c);
            if is_c(k_c,b_c,x_c)==1 && profit_val-b_val+theta*(1-delta)*k_val>=0
                jstar = left_arr(k_c,b_c,x_c); 
                omega = omega_arr(k_c,b_c,x_c);
                kp_c = pol_kp_ind(k_c,b_c,x_c);
                kp_val    = k_grid(kp_c);
                %b_grid_kp = b_grid(kp_c,:)';
                bprime = bprime_arr(k_c,b_c,x_c);
                %v0_int = myinterp1(b_grid_kp,EVx(kp_c,:)',bprime,1);
                v0_int = omega*EVx(kp_c,jstar)+(1-omega)*EVx(kp_c,jstar+1);
                val_c_new(k_c,b_c,x_c) = q*(psi*(theta*(1-delta)*kp_val-bprime)+(1-psi)*v0_int);
            else
                % unconstrained or forced liquidation: val_c is irrelevant
                val_c_new(k_c,b_c,x_c) = theta*(1-delta)*k_val-b_val;
            end
        end %k
    end %b
end %x

% Update
val_c = val_c_new;

end %close howard iterations h_c

end %end function "fun_howard"
%=========================================================================%

function [jl,omega] = myfind_loc(x_grid,xi)

%-------------------------------------------------------------------------%
% DESCRIPTION
% Find jl s.t. x_grid(jl)<=xi<x_grid(jl+1)
% for jl=1,..,N-1
% omega is the weight on x_grid(jl) so that
% omega*x_grid(jl)+(1-omega)*x_grid(jl+1)=xi
% IMPORTANT: omega is not restricted in [0,1]
% INPUTS
% x_grid must be a strictly increasing column vector (nx,1)
% xi must be a scalar
% OUTPUTS
% jl: Left point (scalar)
% omega: weight on the left point (scalar)
% NOTES
% See find_loc_vec.m for a vectorized version.
% See find_loc.m
%-------------------------------------------------------------------------%

nx = size(x_grid,1);

jl = max(min(locate(x_grid,xi),nx-1),1);
%Weight on x_grid(j)
omega = (x_grid(jl+1)-xi)/(x_grid(jl+1)-x_grid(jl));

% IMPORTANT
%omega = max(min(omega,1),0);

end %end function "myfind_loc"

