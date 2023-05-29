function [model_mom_trans,irf] = fun_targets_tran(data_mom_trans,agg_tran,path,agg_ss,prices,calibWeightsTran)
%Purpose: given transition paths, we compute some moments (i.e. drop in
%consumption etc)
%-------------------------------------------------------------------------%
% INPUTS:
% data_mom_trans: Array with data moments transition, (9,4) 
% agg_tran:       Struct with aggregates on transition 
% path:           Struct with aggregates transition (outside loop)
% agg_ss:         Struct with steady-state aggregates 
% prices:         Struct with prices in steady-state
% calibWeightsTran: Array with same dim as data_mom_trans
%
% OUTPUTS:  
% model_mom_trans: Array with model moments transition, same dim as data_mom_trans
% irf:            Struct with IRFs (perc. change rel. to s.s.) on transition
%-------------------------------------------------------------------------%


irf = struct();
% Impulse response function of var in percentage terms
% length is T+1
%---------------- Outer loop ----------------------------%
irf.C_agg    = (path.C-agg_ss.C_agg)/agg_ss.C_agg;              % consumption
irf.KL_ratio = (path.KL_ratio-agg_ss.K_corp/agg_ss.L_corp)/(agg_ss.K_corp/agg_ss.L_corp); % output small firms
irf.q        = (path.q-prices.q)/prices.q;                      % firms discount factor
irf.w        = (path.w-prices.wage)/prices.wage;                % wage

%------------------- Heterog agents --------------------------%
irf.K_agg   = (agg_tran.K_agg-agg_ss.K_agg)/agg_ss.K_agg;       % household capital
irf.K_small = (agg_tran.K_small-agg_ss.K_small)/agg_ss.K_small; % capital small firms
irf.K_corp  = (agg_tran.K_corp-agg_ss.K_corp)/agg_ss.K_corp;    % corporate capital
irf.mass_small = (agg_tran.mass_small-agg_ss.mass_small)/agg_ss.mass_small; % mass of small firms
irf.Y_agg   = (agg_tran.Y_agg-agg_ss.Y_agg)/agg_ss.Y_agg;       % total output
irf.Y_corp  = (agg_tran.Y_corp-agg_ss.Y_corp)/agg_ss.Y_corp;    % corporate output
irf.output_small = (agg_tran.output_small-agg_ss.output_small)/agg_ss.output_small; % output small firms
irf.L_agg   = (agg_tran.L_agg-agg_ss.L_agg)/agg_ss.L_agg;       % total employment
irf.L_small = (agg_tran.L_small-agg_ss.L_small)/agg_ss.L_small; % employment small firms
irf.L_corp  = (agg_tran.L_corp-agg_ss.L_corp)/agg_ss.L_corp;    % corporate employment
irf.entry   = (agg_tran.entry_vec-agg_ss.entry)/agg_ss.entry;   % Measure of entrants
irf.entry_rate   = (agg_tran.entry_rate_vec-agg_ss.entry_rate)/agg_ss.entry_rate;   % entry rate
irf.entry_cost   = (agg_tran.entry_cost_vec-agg_ss.entry_cost)/agg_ss.entry_cost;   % Total cost of entry
irf.exit_rate = (agg_tran.exit_rate_vec-agg_ss.exit_rate)/agg_ss.exit_rate;% exit rate
irf.exit = (agg_tran.exit_vec-agg_ss.exit)/agg_ss.exit;         % Measure of exiting firms
irf.liq     = (agg_tran.liq_vec-agg_ss.liq)/agg_ss.liq;         % liquidation costs
irf.InvK    = (agg_tran.InvK-agg_ss.InvK)/agg_ss.InvK;          % investment of households
irf.InvK_corp = (agg_tran.InvK_corp-agg_ss.InvK_corp)/agg_ss.InvK_corp; % investment of corporate sector
irf.exit_rate_emp = (agg_tran.exit_rate_emp_vec-agg_ss.exit_rate_emp)/agg_ss.exit_rate_emp;% exit rate

n_data = size(data_mom_trans,1);
model_mom_trans = ones(n_data,4);
model_mom_trans(1,1:4) = 100*irf.Y_agg(1:4);        % change in GDP
model_mom_trans(2,1:4) = 100*irf.C_agg(1:4);        % change in consumption
model_mom_trans(3,1:4) = 100*irf.InvK(1:4);     % change in investment
model_mom_trans(4,1:4) = 100*irf.output_small(1:4); % change in small firm revenues
model_mom_trans(5,1:4) = 100*irf.L_agg(1:4);        % change in employment
model_mom_trans(6,1:4) = 100*irf.L_small(1:4);      % change in employment small firms
model_mom_trans(7,1:4) = 100*irf.L_corp(1:4);       % change in employment corporate sector
model_mom_trans(8,1:4) = 100*irf.exit_rate(1:4);    % change in small firm exit rate
% Change in the annual exit rate from 2019 to 2020
exitrate_2019 = 1-(1-agg_ss.exit_rate)^4;
exitrate_2020 = 1-(1-agg_ss.exit_rate)*(1-agg_tran.exit_rate_vec(1))*(1-agg_tran.exit_rate_vec(2))*(1-agg_tran.exit_rate_vec(3));
model_mom_trans(9,1) = 100*(exitrate_2020-exitrate_2019)/exitrate_2019;    %  change in annual firm exit rate
model_mom_trans(10,1:4) = 100*irf.entry_rate(1:4);    % change in small firm entry rate

width = length('Firm exit rate, annual change:')+3;
fprintf("  \n")
disp("--------------------------------------------")
disp("TRANSITION RESULTS")
disp("--------------------------------------------")

% calibWeightsTran has dim: (7,4)

fprintf("%-*s  %-s     %-s     %-s     \n",width,"Description","Data","Model","Weight")
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in GDP q1:",       data_mom_trans(1,1),model_mom_trans(1,1),calibWeightsTran(1,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in GDP q2:",       data_mom_trans(1,2),model_mom_trans(1,2),calibWeightsTran(1,2))
%fprintf("%-*s %-8.4f  %-8.4f  \n", width,"Drop in GDP q3:",               data_mom_trans(1,3),model_mom_trans(1,3))
%fprintf("%-*s %-8.4f  %-8.4f  \n", width,"Drop in GDP q4:",               data_mom_trans(1,4),model_mom_trans(1,4))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f\n", width,"Drop in consumption:",          data_mom_trans(2,1),model_mom_trans(2,1),calibWeightsTran(2,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in investment:",         data_mom_trans(3,1),model_mom_trans(3,1),calibWeightsTran(3,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in small firms output:", data_mom_trans(4,1),model_mom_trans(4,1),calibWeightsTran(4,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment q1:",      data_mom_trans(5,1),model_mom_trans(5,1),calibWeightsTran(5,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment q2:",      data_mom_trans(5,2),model_mom_trans(5,2),calibWeightsTran(5,2))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment small q1:",data_mom_trans(6,1),model_mom_trans(6,1),calibWeightsTran(6,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment small q2:",data_mom_trans(6,2),model_mom_trans(6,2),calibWeightsTran(6,2))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment corp q1:", data_mom_trans(7,1),model_mom_trans(7,1),calibWeightsTran(7,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment corp q2:", data_mom_trans(7,2),model_mom_trans(7,2),calibWeightsTran(7,2))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in small firm exit, q1:",data_mom_trans(8,1),model_mom_trans(8,1),calibWeightsTran(8,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Firm exit rate, annual change:",data_mom_trans(9,1),model_mom_trans(9,1),calibWeightsTran(9,1))
fprintf("%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in small firm entry, q1:",data_mom_trans(10,1),model_mom_trans(10,1),calibWeightsTran(10,1))

fprintf("%-*s  %-8.4f   \n", width,"Total grant amount:",agg_tran.tot_grant);

fprintf("  \n")

end %end function "fun_targets_tran"

