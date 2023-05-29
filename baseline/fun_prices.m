function [prices] = fun_prices(par)
% Purpose: compute prices (q,w,R) in the steady state.
% Usage: [prices] = fun_prices(par)

if isstruct(par)==0
    error('Input par in fun_prices must be a structure!')
end

q        = par.beta;                          % financial discount factor 
FK       = 1/q + par.delta_k -1;              % MPK corporate sector
rental   = FK;                                % rental rate
KL_ratio = fun.optimal_KL(rental,par);        % capital-labor ratio corporate sector
wage     = fun.marg_prod_labor(KL_ratio,par); % real wage

% Pack outputs into a structure:
prices = v2struct(q,rental,wage,KL_ratio);

if (isstruct(prices)==0)
    error('Output <prices> must be a structure')
end

if par.verbose>=1
    fprintf('  \n')
    disp('--------------------------------------------')
    disp('Prices')
    disp('--------------------------------------------')
    fprintf('q:                         %f  \n',q )
    fprintf('K-L Ratio:                 %f  \n',KL_ratio )
    fprintf('Wage:                      %f  \n',wage )
    fprintf(' \n')
end

end %END FUNCTION "fun_prices"
