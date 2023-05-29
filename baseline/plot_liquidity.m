function [ave_b,ave_bk] = plot_liquidity(b_grid,distrib_tran,weights,par)
% DESCRIPTION
% plot_liquidity plots firms average indebtness along the transition
% INPUTS
%   b_grid_nogrant       : Grid for b in transition, (nk,nb,T+1,nn)
%   distrib_tran_nogrant : mu and mu_active, (nk,nb,nx,T+1,nn)
%   weights              : Weights for impact/noimpact and grant/nogrant dim: (nk,nx,nn)
%   par                  : Structure with model parameters
% OUTPUTS
%   ave_b_nogrant        : Avg. b in transition,   dim: (T+1,1)
%   ave_bk_nogrant       : Avg. b/k in transition, dim: (T+1,1)

% Unpack 
nk = par.nk;
nb = par.nb;
nx = par.nx;
nn = par.nn;
T  = par.T;
k_grid = par.k_grid; % (nk,1)
mu_active = distrib_tran.mu_active; %(nk,nb,nx,T+1,nn)

if ~isequal(size(b_grid),[nk,nb,T+1,nn])
    error('b_grid_nogrant not right size')
end

% Compute average b and average b/k in no grant economy, for each t=1,..,T
ave_b  = zeros(T+1,1);
ave_bk = zeros(T+1,1);
for t=1:T+1
    % Initialize for running sum
    ave_b(t)  = 0;
    ave_bk(t) = 0;
    for n_c = 1:nn
    for x_c = 1:nx
        for b_c = 1:nb
            for k_c=1:nk
                k_val = k_grid(k_c);
                b_val = b_grid(k_c,b_c,t,n_c);
                ave_b(t)  = ave_b(t) + b_val*weights(k_c,x_c,n_c)*mu_active(k_c,b_c,x_c,t,n_c);
                ave_bk(t) = ave_bk(t) + (b_val/k_val)*weights(k_c,x_c,n_c)*mu_active(k_c,b_c,x_c,t,n_c);
            end
        end
    end
    end
    % Normalize by mass of active firms in t
    mu_t      = squeeze(sum(mu_active(:,:,:,t,:),2)); % dim:(nk,nx,nn)
    ave_b(t)  = ave_b(t)/sum(weights.*mu_t,'all');
    ave_bk(t) = ave_bk(t)/sum(weights.*mu_t,'all');
end %end t
    
    
end %end function "plot_liquidity"

