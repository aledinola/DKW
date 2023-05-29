function [phi_dist] = gen_phi_dist(bk0_vec,bk0_prob,k_grid,b_grid,x0_prob,prob_k)

% PURPOSE
% gen_phi_dist creates the initial distribution for entrants, Phi(k,b,x)
% INPUTS:
% bk0_vec  : 3*1 vector with values of initial debt-to-asset ratio
% bk0_prob : 3*1 vector with probabilities of initial debt-to-asset ratio
% OUTPUTS
% phi_dist: 3-dim array, initial prob of entrants on (k,b,x)

nk  = size(k_grid,1);
nb  = size(b_grid,2);
nx  = length(x0_prob);
nbk = length(bk0_vec);

if nbk~=length(bk0_prob)
    error('bk0_vec and bk0_prob have different size!')
end
chk = sum(bk0_prob);
if abs(chk-1)>1e-10
    error('bk0_prob does not sum to one!')
end

phi_dist = zeros(nk,nb,nx);

% b0 (scalar) is the initial debt level for entrants
for k_c = 1:nk
    b_gridk = b_grid(k_c,:)';
    kappa = k_grid(k_c);
    for x_c = 1:nx
        for bk_c = 1:nbk
            bk0 = bk0_vec(bk_c);
            b0  = bk0*kappa;
            [b0_ind,omega] = find_loc(b_gridk,b0);
            % unconstraint new entrants
            phi_dist(k_c,b0_ind,x_c)   = omega*bk0_prob(bk_c)*x0_prob(x_c);
            phi_dist(k_c,b0_ind+1,x_c) = (1-omega)*bk0_prob(bk_c)*x0_prob(x_c);
        end
    end
    phi_dist(k_c,:,:) = prob_k(k_c)*phi_dist(k_c,:,:)/sum(phi_dist(k_c,:,:),'all');
end

phi_dist = phi_dist/sum(phi_dist,'all');

validateattributes(phi_dist, {'double'}, {'finite', 'nonnan', 'nonempty','>=',0,'real','size', [nk,nb,nx]})


end %end function "gen_phi_dist"