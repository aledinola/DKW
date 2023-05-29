function [par,bounds_shocks,data_mom_trans,calibWeightsTran] = set_shocks(par,file_shocks)

% In this function we set the shocks parameters, the targets for the
% transition and the calibration weights.

if ~isstruct(par)
    error('input argument "par" must be a structure')
end
if ~ischar(file_shocks)
    error('input argument "file_shocks" must be a character string')
end

%% Read pandemic shock parameters

fprintf('Reading transition shocks from file "%s" \n',file_shocks)
FID = fopen(fullfile(par.InpDir,file_shocks));
C   = textscan(FID,'%s %f');
fclose(FID);
%names  = C{1};
values = C{2};


par.eta_i     = values(1); % fraction of impacted small firms
par.v_corp        = values(2); % TFP shock corporate sector
par.util_shift    = values(3); % Demand shock
par.lsupply_shift = values(4); % Labor supply shock
par.rho_shock     = values(5); % Autocorrelation of shocks

par.v_small       = -1;%% TFP shock on impacted small firms (shutdown shock = -100% impact)
par.v_small_unimp = -0.0; % TFP shock on unimpacted small firms

% Credit shock. Recall that lam=lam0*theta*(1-delta)
par.lambda_shift = -0.158;  % shock to lambda0

% Entry shock
par.mass_shift = 0.24;

% placebo
% par.v_corp = 0;
% par.util_shift = 0;
% par.lsupply_shift = 0;
% par.v_small = 0;
% par.mass_shift = 0;
% par.lambda_shift = 0;


% Bounds for calibrated shocks
bounds_shocks      = zeros(5,2); 

bounds_shocks(1,:) = [0, 1]; % fraction of impacted small firms
bounds_shocks(2,:) = [-0.2, 0]; % v_corp: TFP corporate sector
bounds_shocks(3,:) = [-0.25, 0];  % util_shift: demand shock
bounds_shocks(4,:) = [0, 0.25];   % lsupply_shift: labor supply shock
bounds_shocks(5,:) = [0.0, 0.5];  % rho_shock

% A is TFP in small firms sector or corporate sector
% The second dimension of A_small indicates impactedness:
% i => (1=impacted,2=no impacted)
% margutil is the marginal utility of consumption multiplier. It drops below 1 in the
% first period, then returns gradually to 1.
% lsupply is the marginal utility of leisure multiplier. It increases above 1 in the
% first period, then returns to 1.
par.A_small  = ones(par.T+1,par.ni); % second dim: 1=impacted, 2=unimpacted
par.A_corp   = ones(par.T+1,1);
par.margutil = ones(par.T+1,1);
par.lsupply  = ones(par.T+1,1);
par.lambda_vec = ones(par.T+1,1);
par.mass_vec = ones(par.T+1,1);
% Initialize first period
par.A_small(1,1) = 1+par.v_small; % impacted
par.A_small(1,2) = 1+par.v_small_unimp; % unimpacted

par.A_corp(1)    = 1+par.v_corp;
par.margutil(1)  = 1+par.util_shift;
par.lsupply(1)   = 1+par.lsupply_shift;
par.lambda_vec(1) = 1+par.lambda_shift;
par.mass_vec(1) = 1+par.mass_shift;

for t=2:par.T+1
    par.A_small(t,1) = 1+par.rho_shock^(t-1)*par.v_small; % impacted
    par.A_small(t,2) = 1+par.rho_shock^(t-1)*par.v_small_unimp; % unimpacted

    par.A_corp(t)    = 1+par.rho_shock^(t-1)*par.v_corp;
    par.margutil(t)  = 1+par.rho_shock^(t-1)*par.util_shift;
    par.lsupply(t)   = 1+par.rho_shock^(t-1)*par.lsupply_shift;
    par.lambda_vec(t) = (1+par.rho_shock^(t-1)*par.lambda_shift);
    par.mass_vec(t) = (1+par.rho_shock^(t-1)*par.mass_shift);

end
par.lambda_vec = par.lambda_vec *par.lambda0*par.theta*(1-par.delta_k);
par.mass_vec      = par.mass_vec * par.mass;

%% Load data moments for the transition

% Targets if a period = a quarter.
% 7 variables x 4 quarters: first quarter is 2020Q2
data_mom_trans = nan(10,4);
data_mom_trans(1,:) = [-10.857,-2.246,-0.774,1.256]; % Output of nonfarm business sector
data_mom_trans(2,:) = [-9.667,-1.488,-0.665,2.061]; % Consumption
data_mom_trans(3,:) = [-15.398,-1.723,3.843,3.242]; % Investment
data_mom_trans(4,:) = [-15.650,-4.107,-1.404,0]; % Output of small firms
data_mom_trans(5,:) = [-12.85,-7.578,-5.437,-5.052]; % Employment
% Employment drop by firm size is based on data provided by Cajner et al
% (2020), based on ADP data. Only 2020Q2 is available.
% Ref: sheet "Emp_drop_by_firmsize" of the excel file
% "~/SHARED/Data/ADP/ADP_Cajner_etal_2021.xlsx".
data_mom_trans(6,:) = [-16.02133603,0,0,0]; % Employment small firms
data_mom_trans(7,:) = [-13.24832078,0,0,0]; % Employment large firms
% small firm exit rate change
data_mom_trans(8,:) = [37.84,-0.06,-10.40,-13.85]; % 0 means missing
% change in annual small firm exit rate
data_mom_trans(9,:) = [3.4,0,0,0];% 0 means missing
% small firm entry rate change
data_mom_trans(10,:) = [-12.5,6.25,9.38,12.5]; 
% The following employment drop by firm size comes directly from ADP aggregate data.
% Ref: "~/SHARED/Data/ADP/ADP_NER_History_2021_10.xlsx".
% We cannot control for transitions between firm-size bins using the
% aggregate data. That is, in the pandemic recession
%, large firm may cut
% employment and transition to a smaller firm-size bin. Thus, the
% employment drop in smaller firm-size bins may be underestimated because
% there may be more firms belonging to these bins.
%data_mom_trans(6,:) = [-12.539,-8.306,-6.931,-6.201]; % Employment small firms
%data_mom_trans(7,:) = [-12.068,-8.897,-8.246,-8.131]; % Employment large firms

% Firm exit and entry data are from BLS-BED, see "Facts from BLS.pdf" in
% ~/data/BED/

%% Set calibration weights for the transition
calibWeightsTran      = zeros(10,4);
calibWeightsTran(1,1) = 1.0; % weight on delta GDP q1
calibWeightsTran(1,2) = 1.0; % weight on delta GDP q2
calibWeightsTran(2,1) = 1.0; % weight on delta consumption
calibWeightsTran(3,:) = 0.0; % weight on investment
calibWeightsTran(4,1) = 0.0; % weight on small firm output
calibWeightsTran(5,1) = 1.0; % weight on employment
calibWeightsTran(6,1) = 1.0; % weight on employment small firms
calibWeightsTran(7,1) = 0.0; % weight on employment large firms
calibWeightsTran(8,1) = 0.0; % weight on small firm exit
calibWeightsTran(9,1) = 0.001; % weight on small firm annual exit
calibWeightsTran(10,1) = 1.0; % weight on small firm entry rate

end %end function "set_shocks"