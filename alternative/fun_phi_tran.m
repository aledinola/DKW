function [phi_dist] = fun_phi_tran(phi_dist_ss,b_grid_ss,b_grid,b_tilde,par)
%------------------ DESCRIPTION ------------------------------------------%
% Compute distribution of entrants Phi(b,x,kappa) along the transition and
% for impacted/non impacted firms
% The Phi(:,x,kappa) in the steady-state is defined on b_grid_ss(:,kappa)
% but during the transition the b grid changes: b_grid(:,k_c,t,i_c)
% Therefore, we 

% INPUTS
% phi_dist_ss :: dim(nk,nb,nx)
% b_grid_ss   :: dim(nk,nb)
% b_grid      :: dim(nk,nb,T+1,nn)
% b_tilde     :: dim(nk,T+1,nn)
%-------------------------------------------------------------------------%

% Unpack
T  = par.T;
nx = par.nx;
nb = par.nb;
nk = par.nk;
nn = par.nn; % = ni x ns
lambda  = par.lambda;
k_grid = par.k_grid; % dim nk

phi_dist = zeros(nk,nb,nx,T+1,nn);

% To simplify, we assume that phi_dist along the transition is equal to the
% steady-state phi_dist

% for t = 1:T+1
%     for i_s = 1:ns
%         phi_dist(:,:,:,t,i_s) = phi_dist_ss;
%     end
% end

% For periods 1 to T+1, project Phi in s.s. to the grid in each transition period.
% Note that eta share of firms are impacted.

for n_c = 1:nn
    for t = 1:T+1
        for k_c = 1:nk
            kappa = k_grid(k_c);
            if b_tilde(k_c,t,n_c)>=lambda*kappa
                for x_c = 1:nx
                    phi_dist(k_c,1,x_c,t,n_c) = sum(phi_dist_ss(k_c,:,x_c));
                end
            else
                for b_c = 1:nb
                    if b_grid_ss(k_c,b_c) <= b_tilde(k_c,t,n_c)
                        % If b on the s.s. b_grid is less than b_tilde in the
                        % transition, then put all mass at the s.s. Phi(b,x) on
                        % Phi(1,x,t) of the transition.
                        for x_c = 1:nx
                            phi_dist(k_c,1,x_c,t,n_c)= phi_dist(k_c,1,x_c,t,n_c) + phi_dist_ss(k_c,b_c,x_c);
                        end
                    else
                        % If b on the s.s. b_grid is between b_tilde(t) and theta*kappa,
                        % distribute mass of Phi(b,x) to the bracketing points.
                        b_ss = b_grid_ss(k_c,b_c);
                        [left_loc,omega] = find_loc(squeeze(b_grid(k_c,:,t,n_c))',b_ss);
                        %left_loc = locate_equi(b_grid(:,k_c,t,n_c),b_ss);
                        %left_loc = max(min(locate(b_grid(:,k_c,t,i_c),b_ss),nb-1),1);
                          
                        %Weight on left_loc
                        %omega = (b_grid(left_loc+1,k_c,t,n_c)-b_ss)/(b_grid(left_loc+1,k_c,t,n_c)-b_grid(left_loc,k_c,t,n_c));
                        %omega = max(min(omega,1),0);
                        for x_c = 1:nx
                            % Distribute mass according to the weight omega
                            phi_dist(k_c,left_loc,x_c,t,n_c) =  phi_dist(k_c,left_loc,x_c,t,n_c)+omega*phi_dist_ss(k_c,b_c,x_c);
                            phi_dist(k_c,left_loc+1,x_c,t,n_c) =  phi_dist(k_c,left_loc+1,x_c,t,n_c)+(1-omega)*phi_dist_ss(k_c,b_c,x_c);
                        end
                    end %end if
                end %end b
            end %end if
        end %end kappa
    end %end t
end %end n_c
    

end %END FUNCTION "fun_phi_tran"

