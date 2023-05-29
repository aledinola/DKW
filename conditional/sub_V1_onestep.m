function [V2] = sub_V1_onestep(V1,pi_x,profit_mat,k_grid,q,theta,delta,psi,do_howard)
% DESCRIPTION
% One-step Bellman operator associated to value function V(k,x) for
% unconstrained firms. During the transition, set do_howard=0.
% INPUTS
% "V1"         Value function V(k,x),        dim: (nk,nx)
% "pi_x"       Transition matrix prob(x,x'), dim: (nx,nx)
% "profit_mat" Static profit pi(k,x),        dim: (nk,nx) 
% "k_grid"     Fixed grid for capital,       dim: (nk,1)
% "q,theta,delta,psi" Parameters,            dim: scalars
% OUTPUTS
% "V2"         Value function V(k,x),        dim: (nk,nx)
% NOTES:
% We do maximization over k' on the grid, see code for golden
% that I commented out.

n_howard = 50;

[nk,nx] = size(V1);

V2       = zeros(nk,nx);
kpol_ind = ones(nk,nx);  % Policy function k'(k,x), indexes
kprime   = k_grid;  % Next-period capital: we vectorize over k'. dim:(nk,1)
k_today  = k_grid'; % Current period capital 
%lb = k_grid(1);
%ub = k_grid(nk);

% Compute expected value
V1_max = max(theta*(1-delta)*kprime,V1); % dim: (k,x)
EV = V1_max*pi_x'; % dim: (k,x)

for x_c = 1:nx
    EV_x     = EV(:,x_c);
    profit_x = profit_mat(:,x_c)';
    %EVx = zeros(nk,1);
    %for xp_c = 1:nx
    %    EVx = EVx+pi_x(x_c,xp_c)*max(theta*(1-delta)*kprime,V1(:,xp_c));
    %end
    %for k_c = 1:nk
        %k_val = k_grid(k_c);
        
        %   [x,fval] = golden(f,a,b,tol,P1,P2,...);
        %[~,fval] = golden(@fun_rhs,lb,ub,1e-6);
        %V2(k_c,x_c) = fval;
        RHS = profit_x-fun.adjcost(kprime,k_today,theta,delta)+...
            q*(psi*theta*(1-delta)*kprime+(1-psi)*EV_x);
        [V2(:,x_c),kpol_ind(:,x_c)] = max(RHS,[],1);
        
    %end %k
end % x

%----------------- HOWARD ---------------------------------------------%
% Update V2 a number of times without doing the maximization
if do_howard==1
    for h_c = 1:n_howard
        % Compute expected value
        EVh = max(theta*(1-delta)*kprime,V2)*pi_x';
        for x_c = 1:nx
            EV_x = EVh(:,x_c);
            for k_c = 1:nk
                k_val = k_grid(k_c);
                kopt_ind = kpol_ind(k_c,x_c);
                V2(k_c,x_c) = profit_mat(k_c,x_c)-fun.adjcost_scal(k_grid(kopt_ind),k_val,theta,delta)+...
                    q*(psi*theta*(1-delta)*k_grid(kopt_ind)+(1-psi)*EV_x(kopt_ind));
            end %k
        end % x
    end %end howard h_c
end %end do_howard

% V1_new = ones(nx,nk);
% for k_c = 1:nk
%    kappa = k_grid(k_c);
%    nextV = max(theta*kappa,V1(:,k_c)); %dim: (nx,1)
%    % Continuation value EV := E_x{ max(theta*kappa,V(x') }
%    EV = pi_x*nextV; %dim: (nx,1)
%    % Bellman equation
%    V2(:,k_c) = profit_vec(:,k_c) + q*(psi*theta*kappa+(1-psi)*EV); %dim: (nx,1)
% end

%-------------------- NESTED FUNCTIONS-----------------------------------%

%     function F = fun_rhs(kprime_val)
%
%         EVx_int = myinterp1(k_grid,EVx,kprime_val,1);
%         F = profit_mat(k_c,x_c)-fun.adjcost_scal(kprime_val,k_val,theta)+q*EVx_int;
%
%     end %end nested function

end %end main function

