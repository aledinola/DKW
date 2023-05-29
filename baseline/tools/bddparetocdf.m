function cdf = bddparetocdf(emin, emax, shape, x)

% bddparetocdf calculates the cumulative density (or mass) at x for a
% discretized bounded Pareto distribution over (emin,emax)
% Authors: In Hwan Jo and Tatsuro Senga, later modified by 
% Alessandro Di Nola

paretocdf = 1.0 - (emin/x)^shape;
cdf = paretocdf/(1.0 - (emin/emax)^shape);

end %end function "bddparetocdf"