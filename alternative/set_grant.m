function [par] = set_grant(par)

% DESCRIPTION
% We add new variables to struct "par".
% Set rescue policy parameters Xp and eta for transition
% par.weights has dim: (nk,nx,nn). For each (k,x), the value of weights
% represent the prob. of being in a certain {impact*grant} state.
% nn is a dummy that takes 4 values:
% 1 = grant and impact
% 2 = grant and no impact
% 3 = no grant and impact
% 4 = no grant and no impact
%
% For example, UNTARGETED GRANT (i.e. baseline grant)
% eta_i = prob of being impacted, dim: scalar 
% eta   = prob of receiving the grant, dim: (nk,nx)
% In this case, weights are the same for all (k,x) and since the two events
% (grant and impact) are independent we have:
% weights(:,:,n=1) is eta*eta_i
% weights(:,:,n=2) is eta*(1-eta_i)
% etc.
% NOTE: the targeted grant of the paper is called here "targslim".

if ~isstruct(par)
    error('input argument "par" must be a structure')
end

par.T_grant = 1; % number of periods grant is spread

% baseline amount of grant: 2.5 times the monthly payroll, depending on x in the absence of the shock
if par.grant_flag == 0 % no grant
    par.Xp = 0;
    par.eta = 0.76*ones(par.nk,par.nx); % share of small firms with grant
    
elseif par.grant_flag == 1 % baseline grant
    par.Xp = 2.5/3;
    par.eta = 0.76*ones(par.nk,par.nx); % share of small firms with grant
elseif par.grant_flag == 2 % grant based on x
    par.Xp = 2.5/3;
    par.eta = zeros(par.nk,par.nx);
    for k_c = 1:par.nk
        for x_c = 1:par.nx
            if par.x_grid(x_c)>=1.0 && par.x_grid(x_c)<=2.0
                par.eta(k_c,x_c) = 0.76;
            end
        end
    end
elseif par.grant_flag == 3 % uniform grant large
    par.Xp = 20;%*2.5/3;
    par.eta = 0.76*ones(par.nk,par.nx); % share of small firms with grant
elseif par.grant_flag == 4 % uniform grant small
    par.Xp = 0.5*2.5/3;
    par.eta = 0.76*ones(par.nk,par.nx); % 
elseif par.grant_flag == 5 % only impacted firms get grant
    par.Xp = 2.5/3;
    par.eta = par.eta_i*ones(par.nk,par.nx); % eta==eta_i
elseif par.grant_flag == 6 % only impacted firms get grant, large grant
    par.Xp  = 2.5/3*6.32; % increase Xp 
    par.eta = par.eta_i*ones(par.nk,par.nx); % eta==eta_i
elseif par.grant_flag == 7 % only impacted firms get grant, small grant
    par.Xp  = 2.5/3*0.5; % decrease Xp
    par.eta = par.eta_i*ones(par.nk,par.nx); % eta==eta_i
end


% Whether grant is targeted to impacted firms
par.weights = zeros(par.nk,par.nx,par.nn);
if par.grant_target == 1 % targeted grant
    if max(par.eta,[],"all")<par.eta_i
        disp("Too many impacted firms, not enough grant to target!")
        keyboard
    end
    % targeted to impacted firms
    % i.e. impacted firms receive grant with prob=1 
    %   unimpacted firms receive grant with prob = (eta-eta_i)/(1-eta_i)
    eta_unimp = (par.eta-par.eta_i)./(1-par.eta_i);
    par.weights(:,:,1) = par.eta_i; % impacted, grant
    par.weights(:,:,2) = (1-par.eta_i)*eta_unimp; % unimpacted, grant
    par.weights(:,:,3) = 0; % impacted, nogrant
    par.weights(:,:,4) = (1-par.eta_i) * (1-eta_unimp); % unimpacted, no grant
elseif par.grant_target == 2 % slim target
    if max(par.eta,[],"all")<par.eta_i
        disp("Too many impacted firms, not enough grant to target!")
        keyboard
    end
    % targeted to only impacted firms
    % i.e. impacted firms receive grant with prob = 1 
    %   unimpacted firms receive grant with prob = 0
    %eta_unimp = (par.eta-par.eta_i)./(1-par.eta_i);
    par.weights(:,:,1) = par.eta_i; % impacted, grant
    par.weights(:,:,2) = 0; % unimpacted, grant
    par.weights(:,:,3) = 0; % impacted, nogrant
    par.weights(:,:,4) = (1-par.eta_i); % unimpacted, no grant
else % untargeted (i.e. baseline grant)
    par.weights(:,:,1) = par.eta_i * par.eta;
    par.weights(:,:,2) = (1-par.eta_i)*par.eta;
    par.weights(:,:,3) = par.eta_i * (1-par.eta);
    par.weights(:,:,4) = (1-par.eta_i) * (1-par.eta);
end

% Check: 
% for each x, k, the sum of weights over n_c should be equal to 1
if abs(max(sum(par.weights,3),[],'all')-1)>1e-10 || abs(min(sum(par.weights,3),[],'all')-1)>1e-10
    disp("weights do not sum to 1")
end

end %end function "set_grant"