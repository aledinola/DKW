function [CEV_baseline,CEV_targslim] = fun_welfare(par,C_ss,L_ss,C_baseline,L_baseline,...
    C_nogrant,L_nogrant,C_targslim,L_targslim,margutil,lsupply)

% DESCRIPTION
% We compute the consumption equivalant variation (CEV) of the rescue
% policies relative to the "no grant" economy.
%
% INPUTS
%    C_ss,L_ss:             Consumption and labor in the steady-state 
%    C_baseline,L_baseline: C and L in transition baseline grant, (T+1,1)
%    C_nogrant,L_nogrant:   C and L in transition no grant, (T+1,1)
%    C_targslim,L_targslim: C and L in transition targeted grant, (T+1,1)
%    margutil,lsupply:      paths for demand and labor supply shocks, (T+1,1)
% OUTPUTS
%    CEV_baseline: CEV of baseline grant vs no grant
%    CEV_targslim: CEV of targeted grant vs no grant

if ~isstruct(par)
    error('Input argument par must be a struct')
end
if ~isequal(size(C_baseline),[par.T+1,1])
    error('C_baseline must be a column vector with T+1 elements')
end

% Compute value function for the representative household
V_baseline     = 0;
V_nogrant      = 0;
V_targslim     = 0;
V_nogrant_cons = 0;

for t = 1:par.T+1
    V_baseline = V_baseline + par.beta^(t-1)* ...
        fun.utility(C_baseline(t),margutil(t),L_baseline(t),lsupply(t),par);
    V_nogrant = V_nogrant + par.beta^(t-1)* ...
        fun.utility(C_nogrant(t),margutil(t),L_nogrant(t),lsupply(t),par);
    V_targslim = V_targslim + par.beta^(t-1)* ...
        fun.utility(C_targslim(t),margutil(t),L_targslim(t),lsupply(t),par);
    % Here we consider only the consumption part of the value function in
    % the first four quarters
    if t<=4
        V_nogrant_cons = V_nogrant_cons + par.beta^(t-1)* ...
            fun.utility_consumption(C_nogrant(t),margutil(t),par);
    end
end
% Compute continuation value after T+1 periods. Recall that in the s.s. the
% shocks are equal to one.
V_ss      = 1/(1-par.beta) * fun.utility(C_ss,1,L_ss,1,par);

% Add continuation value
V_baseline     = V_baseline + par.beta^(par.T+1)*V_ss;
V_nogrant      = V_nogrant + par.beta^(par.T+1)*V_ss;
V_targslim     = V_targslim + par.beta^(par.T+1)*V_ss;
%V_nogrant_cons = V_nogrant_cons + par.beta^(par.T+1)*V_ss_cons;

% Compute CEV
CEV_baseline = ((V_baseline - V_nogrant)/V_nogrant_cons+1)^(1/(1-par.sigma))-1;
CEV_targslim = ((V_targslim - V_nogrant)/V_nogrant_cons+1)^(1/(1-par.sigma))-1;

end %end function "fun_welfare"