function [mu1] = sub_mu_onestep(mu,phi_dist,pol_kp_ind,pol_exit,pol_entry,...
    left_loc_arr,omega_arr,pi_x,mass,psi)

% DESCRIPTION:
% It does one step of the operator over the distribution mu^0
% This function is called either in the steady-state by fun_distrib*. For
% speed reasons we precomputed the indexes and weights for interpolation in
% left_loc_arr and omega_arr.
% INPUTS:
%   "mu"            Initial distribution mu^0(k,b,x)
%   "phi_dist"      Distribution of entrants Phi(k,b,x)
%   "pol_kp_ind"    Policy k'(k,b,x), indexes
%   "pol_exit"      Policy for exit d^l(k,b,x)
%   "pol_entry"     Policy for exit d^e(k,b,x)
%   "left_loc_arr"  Indexes for interpolation of b'(k,b,x)
%   "omega_arr"     Weights for interpolation of b'(k,b,x)
%
% OUTPUT
%   "mu1"         Updated distribution mu^0(k,b,x)

[nk,nb,nx]  = size(mu);

mu1 = zeros(nk,nb,nx);


for x_c = 1:nx % current productivity
    for b_c = 1:nb % current debt
        for k_c = 1:nk % current capital
            dexit = psi+(1-psi)*pol_exit(k_c,b_c,x_c); %1=exit,0=stay
            entry = pol_entry(k_c,b_c,x_c);
            knext_ind = pol_kp_ind(k_c,b_c,x_c);
            left_loc = left_loc_arr(k_c,b_c,x_c);
            %Weight on left_loc
            omega = omega_arr(k_c,b_c,x_c);
            %for xp_c = 1:nx
            mu1(knext_ind,left_loc,x_c)   = mu1(knext_ind,left_loc,x_c)+omega*(1-dexit)*mu(k_c,b_c,x_c)+ ...
                omega*mass*entry*phi_dist(k_c,b_c,x_c);
            mu1(knext_ind,left_loc+1,x_c) = mu1(knext_ind,left_loc+1,x_c)+(1-omega)*(1-dexit)*mu(k_c,b_c,x_c)+ ...
                (1-omega)*mass*entry*phi_dist(k_c,b_c,x_c);
            %end %xp_c
        end %k_c
    end %b_c
end %x_c

for k_c = 1:nk
    % Matrix multiplication: mu1(k',b',x)*pi(x,x')==> mu1(k',b',x')
    temp = squeeze(mu1(k_c,:,:));
    mu1(k_c,:,:) = temp*pi_x;
end

end %END FUNCTION <sub_mu_onestep>
