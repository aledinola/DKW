function F = dyn_eqn_capital(k_next,k_current,par,t)
% Difference equation for capital to solve during the transition

k_next = max(k_next,1e-10);

lhs = par.beta*k_current^par.alpha;
%rhs1 = par.A_corp(t+1)*par.margutil(t)/par.A_corp(t);
rhs1 = par.A_corp(t+1)*par.margutil(t)*par.lsupply(t)/(par.A_corp(t)*par.lsupply(t+1));
rhs2 = k_next^par.alpha/(1-par.delta_k+par.alpha*par.A_corp(t+1)*par.A*k_next^(par.alpha-1));


F = lhs-rhs1*rhs2;


end