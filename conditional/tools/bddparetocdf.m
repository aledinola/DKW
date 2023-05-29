function cdf = bddparetocdf(emin, emax, shape, x)

% This calculates the cumulative density (or mass) at x for a
% discretized bounded Pareto distribution over (emin,emax)

paretocdf = 1.0 - (emin/x)^shape;
cdf = paretocdf/(1.0 - (emin/emax)^shape);

end %end function "bddparetocdf"