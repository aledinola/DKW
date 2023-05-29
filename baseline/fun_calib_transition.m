function [distance] = fun_calib_transition(x_in,par,data_mom_trans,calibWeightsTran,bounds_shocks)

global obj_tran_best


%check if parameters are within lower bounds
if any(x_in<bounds_shocks(:,1))
    warning('Lower bounds on parameters violated!')
    distance = realmax; %1e+10
    return
end
%check if parameters are within upper bounds
if any(x_in>bounds_shocks(:,2))
    warning('Upper bounds on parameters violated!')
    distance = realmax;
    return
end

par.eta_i       = x_in(1);
v_corp        = x_in(2);
util_shift    = x_in(3);
lsupply_shift = x_in(4);
rho_shock     = x_in(5);


% Store shocks into the "par" structure to be passed to fun_transition
par.A_small  = ones(par.T+1,par.ni);
par.A_corp   = ones(par.T+1,1);
par.margutil = ones(par.T+1,1);
par.lsupply  = ones(par.T+1,1);

% Initialize first period
par.A_small(1,1) = 1+par.v_small; % impacted
par.A_corp(1)    = 1+v_corp;
par.margutil(1)  = 1+util_shift;
par.lsupply(1)   = 1+lsupply_shift;

for t=2:par.T+1
    par.A_small(t,1) = 1+rho_shock^(t-1)*par.v_small; % impacted
    par.A_corp(t)    = 1+rho_shock^(t-1)*v_corp;
    par.margutil(t)  = 1+rho_shock^(t-1)*util_shift;
    par.lsupply(t)   = 1+rho_shock^(t-1)*lsupply_shift;
end


% Whether grant is targeted to impacted firms
par.weights = zeros(par.nx,par.nk,par.nn);
if par.grant_target == 1
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

else % untargeted
    par.weights(:,:,1) = par.eta_i * par.eta;
    par.weights(:,:,2) = (1-par.eta_i)*par.eta;
    par.weights(:,:,3) = par.eta_i * (1-par.eta);
    par.weights(:,:,4) = (1-par.eta_i) * (1-par.eta);
end

if par.grant_flag==0
    error('Must do calibration with grant ON')
end

if par.grant_target==1
    error('Must do calibration with UNTARGETED grant')
end

fprintf('============================================================= \n')
fprintf('Start steady-state computation... \n')
fprintf(' \n')
tic
% Exit rate in the s.s.: agg.exit_rate
[sol,agg,b_grid,distribS,prices,model_mom,flag_ss,par] = fun_steady_state(par);
if flag_ss<0
    warning("steady-state not found!")
else
    fprintf('Steady-state computed! \n')
end
toc

fprintf(' \n')
fprintf('Start transition... \n')
fprintf('Transition periods T = %d \n',par.T)
fprintf('grant_flag           = %d \n',par.grant_flag)

tic
% Exit rate in the transition: agg_tran.exit_rate_vec
[agg_tran,path,conv_flag,pol_tran,distrib_tran] = fun_transition(par,sol,agg,distribS,prices,b_grid);
disp('Time to compute transition:')
toc
if conv_flag<0
    warning("Be careful!")
    fprintf('Transition did not converge \n' )
end

% Display shocks
%width = length('Drop in small firms output:');
fprintf("  \n")
disp("--------------------------------------------")
disp("TRANSITION SHOCKS")
disp("--------------------------------------------")
fprintf("eta_i:       %8.6f  \n",par.eta_i)
fprintf("v_small:       %8.6f  \n",v_small)
fprintf("v_corp:        %8.6f  \n",v_corp)
fprintf("util_shift:    %8.6f  \n",util_shift)
fprintf("lsupply_shift: %8.6f  \n",lsupply_shift)
fprintf("rho_shock      %8.6f  \n",rho_shock)

fprintf("  \n")

% Compute moments for transition
%if conv_flag>=0
[model_mom_trans,irf] = fun_targets_tran(data_mom_trans,agg_tran,path,agg,prices,calibWeightsTran);
%end

if numel(data_mom_trans)~=numel(model_mom_trans)
    error('data_mom_trans and model_mom_trans do not conform')
end

if numel(data_mom_trans)~=numel(calibWeightsTran)
    error('data_mom_trans and calibWeightsTran do not conform')
end

n_data = size(data_mom_trans,1);
dev = zeros(n_data,1);
for i=1:n_data
    dev(i) = sum(calibWeightsTran(i,:).*((model_mom_trans(i,:)-data_mom_trans(i,:)) ...
        ./(abs(data_mom_trans(i,:))+1e-20)).^2);
end
distance = sum(dev);

if isnan(distance)
    distance = realmax;
end

fprintf('\n');
fprintf("Sqr. dist GDP:           %f  \n",dev(1) )
fprintf("Sqr. dist C:             %f  \n",dev(2) )
fprintf("Sqr. dist Inv:           %f  \n",dev(3) )
fprintf("Sqr. dist Y small firms: %f  \n",dev(4) )
fprintf("Sqr. dist Emp:           %f  \n",dev(5) )
fprintf("Sqr. dist Emp small:     %f  \n",dev(6) )
fprintf("Sqr. dist Emp corp:      %f  \n",dev(7) )
fprintf("Sqr. dist Exit rate:     %f  \n",dev(8) )
fprintf("Sqr. dist Exit rate (annual):     %f  \n",dev(9) )
fprintf("Sqr. dist Entry rate:     %f  \n",dev(10) )

fprintf("Total distance:          %f  \n",distance )
fprintf('=============================================================\n');


if distance<obj_tran_best
    obj_tran_best = distance;
    % append results to txt file
    append_tran_txt(distance,x_in,data_mom_trans,model_mom_trans,calibWeightsTran,par.InpDir);
end

end %END function <fun_calib_transition>

