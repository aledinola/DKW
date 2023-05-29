function [out_zombie] = fun_zombie(pol_tran,distrib_tran,path,par)

% DESCRIPTION
%   Compute statistics of zombie firms (saved and zombie) and non-zombie
%   saved firms. 
% INPUTS
%   pol_tran : fields used: pol_exit(k,b,x,t,n) and V1(k,x,t,n) 
%   distrib_tran : Struct with distrib in transition. Field used:
%   mu(k,b,x,t,n).
%   path     : Struct with prices in the transition. Field used: w(t)
%   par      : Struct with parameters. Fields used include weights(k,x,n).
% OUTPUT
%   out_zombie    : Struct with fields ave_x_zombie, ave_x_no_zombie, ave_b_zombie, ave_b_no_zombie
% ave_k_zombie, ave_k_no_zombie, ave_l_zombie, ave_l_no_zombie, ave_y_zombie, ave_y_no_zombie, mass_zombie
% mass_no_zombie, mass_not_saved, mass_nogrant_exit.


% NOTES
%   Called by main_plots.m


% Create indicator zombie, dim: nb,nx,nk,ni
%   0 = firm is not saved (either always in or always out)
%   1 = firm is saved, and zombie
%   2 = firm is saved, not a zombie

% 1 = impact, grant,  2 = no impact, grant, 3 = impact, no grant, 4 =  no impact, no grant
zombie = zeros(par.nk,par.nb,par.nx,par.ni);
t_c = 1; % consider only the impact period
for i_c = 1:par.ni  % impact vs no impact
    for k_c = 1:par.nk
        for x_c = 1:par.nx
            for b_c = 1:par.nb
                % Check whether the firm is saved
                if pol_tran.pol_exit(k_c,b_c,x_c,t_c,i_c)<0.5 && ...
                    pol_tran.pol_exit(k_c,b_c,x_c,t_c,i_c+2)>=0.5   
                    % if firm does not exit under grant (n_c = i_c) but
                    % exits under no grant (n_c = i_c+2), the firm is saved
                    zombie(k_c,b_c,x_c,i_c) = 2; 
                    % Check whether the firm is zombie
                    % Note that V1 has dimension nx,nk,T+1,nn

                    if (pol_tran.V1(k_c,x_c,t_c,i_c+2)<par.theta*(1-par.delta_k)*par.k_grid(k_c))
                        % if firm's unconstrained value is smaller than
                        % liquidation value, the firm is a zombie
                        zombie(k_c,b_c,x_c,i_c) = 1;
                    end
                end
            end
        end
    end
end

% What fraction of saved firms are zombie? What are the ave. x, k, b of
% zombies and non-zombies?
% Use distribution prior to exit, called mu^0 in the draft
t_c = 1; % consider only the impact period

mass_not_saved = 0;
mass_zombie    = 0;
mass_no_zombie = 0;
mass_nogrant_exit = 0; % firms with grant that would have exited without a grant
mass_allfirms = 0; % mass of all active firms
ave_x_zombie = 0;
ave_b_zombie = 0;
ave_k_zombie = 0;
ave_l_zombie = 0;
ave_y_zombie = 0;
ave_x_no_zombie = 0;
ave_b_no_zombie = 0;
ave_k_no_zombie = 0;
ave_l_no_zombie = 0;
ave_y_no_zombie = 0;
for i_c = 1:par.ni  % impact vs no impact with grant
    for k_c = 1:par.nk
        for x_c = 1:par.nx
            for b_c = 1:par.nb
                % Mass of firms
                mass_zombie = mass_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==1);
                mass_no_zombie = mass_no_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==2);
                mass_not_saved = mass_not_saved+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==0);
                % Mass of firms (with and without the grant) that would have exited without the grant
                mass_nogrant_exit = mass_nogrant_exit+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)  ...
                    * (pol_tran.pol_exit(k_c,b_c,x_c,t_c,i_c+2)>=0.5);
                % Mass of all firms
                mass_allfirms = mass_allfirms + par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c);
                % average productivity
                ave_x_zombie = ave_x_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==1)*par.x_grid(x_c);
                ave_x_no_zombie = ave_x_no_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==2)*par.x_grid(x_c);
                % average debt
                ave_b_zombie = ave_b_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==1)*pol_tran.b_grid(k_c,b_c,t_c,i_c);
                ave_b_no_zombie = ave_b_no_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==2)*pol_tran.b_grid(k_c,b_c,t_c,i_c);
                % average capital
                ave_k_zombie = ave_k_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==1)*par.k_grid(k_c);
                ave_k_no_zombie = ave_k_no_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==2)*par.k_grid(k_c);
                % average labor
                % fun_l(x,wage,kappa,par)
                xval = par.A_small(t_c,i_c)*par.x_grid(x_c);
                wage = path.w(t_c);
                lval = fun.fun_l(xval,wage,par.k_grid(k_c),par);
                ave_l_zombie = ave_l_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==1)*lval;
                ave_l_no_zombie = ave_l_no_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==2)*lval;
                % average output
                % prod_small(x,kappa,l_small,c,par)
                cf = 0;%par.fixcost(k_c);
                yval = fun.prod_small(xval,par.k_grid(k_c),lval,cf,par);
                ave_y_zombie = ave_y_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==1)*yval;
                ave_y_no_zombie = ave_y_no_zombie+par.weights(k_c,x_c,i_c)*distrib_tran.mu(k_c,b_c,x_c,t_c,i_c)*(zombie(k_c,b_c,x_c,i_c)==2)*yval;

            end
        end
    end
end


% normalize by mass of firms
out_zombie.ave_x_zombie = ave_x_zombie / mass_zombie;
out_zombie.ave_x_no_zombie = ave_x_no_zombie / mass_no_zombie;
out_zombie.ave_b_zombie = ave_b_zombie / mass_zombie;
out_zombie.ave_b_no_zombie = ave_b_no_zombie / mass_no_zombie;
out_zombie.ave_k_zombie = ave_k_zombie / mass_zombie;
out_zombie.ave_k_no_zombie = ave_k_no_zombie / mass_no_zombie;
out_zombie.ave_l_zombie = ave_l_zombie / mass_zombie;
out_zombie.ave_l_no_zombie = ave_l_no_zombie / mass_no_zombie;
out_zombie.ave_y_zombie = ave_y_zombie / mass_zombie;
out_zombie.ave_y_no_zombie = ave_y_no_zombie / mass_no_zombie;

out_zombie.mass_zombie = mass_zombie;
out_zombie.mass_no_zombie = mass_no_zombie;
out_zombie.mass_not_saved = mass_not_saved;
out_zombie.mass_nogrant_exit = mass_nogrant_exit;
out_zombie.mass_allfirms = mass_allfirms;

end %end function