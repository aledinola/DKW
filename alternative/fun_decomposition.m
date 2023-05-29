function [Delta] = fun_decomposition(T_last,distribS,prices,...
    distrib_tran,path,weights,par)

% DESCRIPTION
%   Decompose change in output small firmts relative to the steady state
%   deltaY(x) = deltaTFP(x)+deltaL(x)+deltaExit(x)
% INPUTS
%   T_last   : no. of periods we consider
%   distribS : Struct with distrib in steady-state, mu and mu_active(k,b,x)
%   prices   : Struct with prices in s.s., wage
%   distrib_tran : Struct with distrib in transition, mu
%                  mu_active(k,b,x,t,n)
%   path     : Struct with prices in the transition
%   weights  :
%   par      : Struct with parameters
% OUTPUT
%   Delta    : Struct with fields deltaY,deltaTFP,deltaL,deltaExit.
% NOTES
%   Called by main_plots.m

nb      = par.nb;     % no. of grid points for debt, "b"
nx      = par.nx;     % no. grid points for x
nk      = par.nk;     % no. grid points for kappa
ni      = par.ni;     % impacted/unimpacted
ns      = par.ns;     % grant/no grant
nn      = par.nn;
x_grid  = par.x_grid; % grid for productivity , "x"
k_grid  = par.k_grid; % grid for capital small firms, "kappa"
fixcost = par.fixcost; % vector of fixed costs, one for each kappa

% SS dim: (k,b,x)
% Tran dim: (k,b,x,t,impact x grant)
deltaY    = zeros(nk,nx);
deltaTFP  = zeros(nk,nx);
deltaL    = zeros(nk,nx);
deltaExit = zeros(nk,nx); %change in Y due to change in mu (i.e. entry and exit)

for n_c = 1:nn
    [i_c,~] = ind2sub([ni,ns],n_c); % i_c = impact indicator; s_c = grant indicator

    for t_c = 1:T_last
        for x_c = 1:nx
            for b_c = 1:nb % current debt
                for k_c = 1:nk
                    kappa = k_grid(k_c);
                    fixval = fixcost(k_c);
                    x1  = par.A_small(t_c,i_c)*x_grid(x_c);
                    xss = x_grid(x_c); % steady-state
                    %         fun_l(x,wage,kappa,par)
                    l1  = fun.fun_l(x1,path.w(t_c),kappa,par);
                    lss = fun.fun_l(xss,prices.wage,kappa,par);
                    %prod_small(x,kappa,l_small,c,par)
                    y1 = weights(k_c,x_c,n_c)*distrib_tran.mu_active(k_c,b_c,x_c,t_c,n_c)*fun.prod_small(x1,kappa,l1,fixval,par);
                    yss = weights(k_c,x_c,n_c)*distribS.mu_active(k_c,b_c,x_c)*fun.prod_small(xss,kappa,lss,fixval,par);
                    deltaY(k_c,x_c) = deltaY(k_c,x_c)+(y1-yss);
                    deltaTFP(k_c,x_c) = deltaTFP(k_c,x_c)+weights(k_c,x_c,n_c)*distribS.mu_active(k_c,b_c,x_c)*(fun.prod_small(x1,kappa,lss,fixval,par)-fun.prod_small(xss,kappa,lss,fixval,par));
                    deltaL(k_c,x_c) = deltaL(k_c,x_c)+weights(k_c,x_c,n_c)*distribS.mu_active(k_c,b_c,x_c)*(fun.prod_small(x1,kappa,l1,fixval,par)-fun.prod_small(x1,kappa,lss,fixval,par));
                    deltaMU = weights(k_c,x_c,n_c)*distrib_tran.mu_active(k_c,b_c,x_c,t_c,n_c)-weights(k_c,x_c,n_c)*distribS.mu_active(k_c,b_c,x_c);
                    deltaExit(k_c,x_c) = deltaExit(k_c,x_c)+deltaMU*fun.prod_small(x1,kappa,l1,fixval,par);
                end % k
            end % b
        end % x
    end % t_c
end %n_c
check = deltaY-(deltaTFP+deltaL+deltaExit);
if max(abs(check),[],'all')>1e-12
    disp(check)
    keyboard
end

Delta.deltaY    = deltaY;
Delta.deltaTFP  = deltaTFP;
Delta.deltaL    = deltaL;
Delta.deltaExit = deltaExit;

end %end function "fun_decomposition"

