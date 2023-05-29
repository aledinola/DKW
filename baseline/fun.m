classdef fun
    %This class contains all functions related to the model (functional
    %forms)
    
    methods (Static)
        
        function adj = adjcost_scal(kprime,k,theta,delta)
            % Compute capital adjustment costs (see Section 3.1 draft)
            % INPUTS
            % kprime: Next-period capital, scalar
            % k: Current period capital, scalar
            % theta: resale value, in (0,1)
            if kprime>=(1-delta)*k
                adj = kprime-(1-delta)*k;
            else
                adj = theta*(kprime-(1-delta)*k);
            end
        end
        
        function adj = adjcost(kprime,k,theta,delta)
            % Compute capital adjustment costs (see Section 3.1 draft)
            % This function is the vectorized version of adjcost_scal
            % INPUTS
            % kprime: Next-period capital, dim (nk,1)
            % k: Current period capital, scalar or dim (1,nk)
            % theta: resale value, in (0,1)
             adj = kprime-(1-delta)*k;
             adj(kprime<(1-delta)*k) = theta*adj(kprime<(1-delta)*k);
        end

        function c = fun_fixcost(k,fixcost1,fixcost2,par)
            % Compute fixed cost of operation
            % INPUTS
            % k: Current period capital x, scalar
            % fixcost1,fixcost2: parameters of the fixed cost function
            c = fixcost1 + fixcost2 * k;
            %c = fixcost1 + fun.prod_small(1,k,1,0,par)*fixcost2;
        end
        
        function C = C_foc_labor(wage,par)
            % Consumption as a function of the wage
            % from the FOC labor supply household
            
            C = (wage/par.zeta)^(1/par.sigma);
            
        end
        
        function F = prod_corp(KL_ratio,L,par)
            % Production function corporate sector
            
            F = par.A*KL_ratio^par.alpha*L;
            
        end
        
        function F = optimal_KL(rental,par)
            % Optimal KL ratio given rental rate
            % from FOC "rental = marginal product of capital"
            
            F = (rental/(par.A*par.alpha))^(1/(par.alpha-1));
            
        end
        
        function F = KL_tran(wage,A_corp_t,par)
            % Optimal KL ratio from FOC wrt labor demand in corp sector
            
            F = (wage/((1-par.alpha)*A_corp_t*par.A))^(1/par.alpha);
            
        end
        
        function F = marg_prod_labor(KL_ratio,par)
            % Marginal product of labor in corporate sector
            
            F = par.A*(1-par.alpha)*(KL_ratio)^(par.alpha);
            
        end
        
        function F = marg_prod_capital(KL_ratio,par)
            % Marginal product of capital in corporate sector
            
            F = par.A*par.alpha*KL_ratio^(par.alpha-1);
            
        end
        
        function F = prod_small(x,kappa,labor,c,par)
            % Production function small firms sector, depends on
            % productivity "x", capital "kappa", labor and fixed cost c.
            % Note: c is fixed cost that depends on kappa
            A = par.A;
            gamma1=par.gamma1;
            gamma2=par.gamma2;
            F = A*x.*(kappa.^gamma1.*labor.^(1-gamma1)).^gamma2-c;
            
        end
        
        function F = fun_l(x,wage,k,par)
            % Labor demand by small business with productivity x and
            % capital k
            % INPUTS
            % "x" Idios. productivity
            % "k" Capital
            gamma1 = par.gamma1;
            gamma2 = par.gamma2;
            aux = (1-gamma1)*gamma2;
            F = (wage./(par.A*x*aux)).^(1/(aux-1)).*k.^(-gamma1*gamma2/(aux-1));
        end
        
        function F = fun_profit(x,k,c,wage,par)
            % Static profit by small business with productivity x
            % and capital k. The fixed cost c^f is already included in
            % fun.prod_small.
            
            l_small = fun.fun_l(x,wage,k,par);
            F = fun.prod_small(x,k,l_small,c,par) - wage*l_small;
            
        end
        
        function F = dyn_eqn_capital(k_next,k_current,par,t)
            % Difference equation for capital to solve during the transition
            k_next = max(k_next,1e-10);
            
            lhs = par.beta*k_current^par.alpha;
            rhs1 = par.A_corp(t+1)*par.margutil(t)*par.lsupply(t)/(par.A_corp(t)*par.margutil(t+1)*par.lsupply(t+1));
            rhs2 = k_next^par.alpha/(1-par.delta_k+par.alpha*par.A_corp(t+1)*par.A*k_next^(par.alpha-1));
            
            F = lhs-rhs1*rhs2;
            
        end

        function U = utility_consumption(C,D,par)
            % Utility from consumption given consumption C and util shifter
            % D
            U = D.*C.^(1-par.sigma)/(1-par.sigma);

        end

        function U = utility_leisure(L,zeta_shift,par)
            % Utility from leisure given labor supply L and zeta shifter
            % zeta_shift
            U = par.zeta*zeta_shift.*(1-L);

        end

        function U = utility(C,D,L,zeta_shifter,par)
            % Utility given consumption C, util shifter D, labor supply L,
            % and zeta shifter zeta_shifter.
            U = fun.utility_consumption(C,D,par) + fun.utility_leisure(L,zeta_shifter,par);
        end

    end %END METHODS
end %END CLASS <fun>