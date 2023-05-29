%% Di Nola, Kaas, Wang (2021)
%% #VC# V23b

%{
----------------------- LEGEND --------------------------------------------

Rescue Policies for Small Businesses During the Covid-19 Recession
May 2021.
This main script calls <fun_steady_state> to compute the pre-pandemic
steady-state of the model. Then it calls <fun_transition> which computes
the transition of the economy after a one-period pandemic shock hits the
economy in t=1. Recall that t=1 in the code corresponds to t=0 in the
draft.
Added PPP policy in transition
To make plots: "main_plots", "make_plots_compare"
To make tables: "main_tables"

----------------------- END LEGEND ----------------------------------------
%}

clear
clc
close all
format long g
digits(100)
% Add path to numerical tools
addpath(genpath(fullfile('..','tools')));
format long g

disp("Hello V24!")

%% Set flags

par.do_calib = 2; % 0 = steady-state,
% 1 = steady-state+transition,
% 2 = calibration steady-state,
% 3 = calibration of transition shocks

par.FigDir   = fullfile('figures','grant'); % folder to save figures
par.TabDir   = fullfile('tables'); % folder to save tables
par.InpDir   = fullfile('inputs'); % folder to read parameters
do_plots     = 1; %flag 0/1 to draw plots
par.do_plots_ss = 0; % flag 0/1 to draw plots for steady-state
do_save      = 0; %flag 0/1 to save plots as png and mat files
par.do_table = 1; %flag 0/1 to write tables on screen
par.do_tex   = 0; %flag 0/1 to write tables on tex (only if do_table=1)
par.verbose  = 0; %flag 0/1/2: 0 no display at all, 1=moderate, 2=disp everything
do_debug     = 0;
est_algo     = 'fminsearch';%'simulan';%
makeCompleteLatexDocument = 1; %flag 0/1 to generate a stand-alone tex doc
load_results = 1;
file_params  = 'estim_params.txt';
file_excel   = 'calibrationV1.xlsx' ; % store here the results of comp stat

% Check if folders exist and if not create them
if ~(isfolder(par.FigDir))
    disp("Folder does not exist, creating it now..")
    mkdir(par.FigDir)
end
if ~(isfolder(par.TabDir))
    disp("Folder does not exist, creating it now..")
    mkdir(par.TabDir)
end

%% Numerical parameters

par.T            = 240; % length transition
par.nx           = 120;%30; %100; % grid size for productivity
par.nb           = 120;%300;%200; %100; % grid size for debt
par.ns           = 2; %  impacted vs unimpacted
par.tol_vfi      = 1e-12; % tolerance for VFI
par.max_iter     = 4000; % max iter for VFI
par.tol_dist     = 1e-10; % tolerance for distribution
par.maxiter_dist = 4000; % max no iter for distribution
par.max_iter_tr  = 100; % max no iter for transition
par.tol_tran     = 0.0001; % tol for transition
par.dampening    = 0.9; % dampen update in bisection
par.matrixInv    = 0; %if 1, use matrix inversion to solve for steady-state mu

%% Set exogenous parameters (not part of internal calibration)

par.beta     = 0.989; % discount factor
par.alpha    = 0.3;   % Cobb-Douglas exponent on capital
par.delta_k  = 0.015; % depreciation rate capital
par.delta_ks = 0.015; % depreciation rate capital owned by small firms
par.gamma    = 0.6;   % span of control parameter in f(l) small firms
par.A        = 0.25;  % production function shifter for the corporate sector
par.psi      = 0;     % exogenous exit rate
par.fixp     = 0;     % share of rental capital of small firms
%par.theta_aux =  0.5; % liquidation value of small firms. theta = theta_aux*(1-fixp)

%% Initial conditions for internal parameters
% If load_results=1, then read parameters from 'estim_params.txt'
if load_results==0
    par.mass    = 1.64548014305575;     % mass of potential entrants
    % it is \xi in the draft
    par.kappa   =  37.3601857658672 ;   % fixed amount of capital for small firms
    par.zeta    = 2.71012295970552 ;    % utility of leaisure
    
    par.epsx    = 0.0864797673312433 ;  % standard deviation of innovation in AR(1) for x'|x
    par.rhox    = 0.789682605114831;    % persistence in AR(1) for x'|x
    par.x0      = 2.89185159031884;     % Ln(x0) mean of AR(1)
    
    par.lambda = 0.0691451225103525;    % exponential distribution parameter for (theta-b) of potential entrants
    par.xi     = 0.533646822820499;     % productivity gap of potential entrants relative to incumbents
    par.fixp     = 0.0;  % share of rental capital of small firms. 0<fixp<1
    par.fixcost  = 0.0;    % fixed operating cost small firm
    
elseif load_results==1
    fprintf('Reading parameters from file "%s" \n',file_params)
    FID = fopen(fullfile(par.InpDir,file_params));
    C = textscan(FID,'%s %f');
    fclose(FID);
    
    names  = C{1};
    values = C{2};
    
    for i=1:numel(names)
        par.(names{i}) = values(i);
    end
    
end

% Set bounds for internal parameters
% [First number is lower bound, second number is the upper bound]

bounds.mass    = [0.0001, 1000000];       % mass of potential entrants
bounds.theta_aux  = [0.001, 1]; % liquidation value of small firms.
%bounds.fixp    = [0, 0.999];    % share of rental capital of small firms. 0<fixp<1
% it is \xi in the draft
bounds.kappa   = [0.5 10000];
%bounds.psi     = [0.0, 0.05]; % exogenous exit probability
bounds.zeta    = [0.1,10];     % utility of leaisure

bounds.epsx    = [0.01,0.5];   % standard deviation of innovation in AR(1) for x'|x
bounds.rhox    = [0.1, 1]; % persistence in AR(1) for x'|x
bounds.x0      = [1, 5.0];  % Ln(x0) mean of AR(1)

bounds.lambda = [0.001, 5];    % exponential distribution parameter for (theta-b) of potential entrants
bounds.xi     = [0.1, 1.0];    % productivity gap of potential entrants relative to incumbents
bounds.fixcost     = [0.000, 100];    % additional fixed operation cost

% Set parameter names IN THE SAME ORDER AS THEY APPEAR in bounds and guess
% calibNames must be column cell array of characters
calibNames = {'mass';
    'theta_aux';
    %'fixp';
    'fixcost';
    'kappa';
    'zeta';
    %'psi';
    'epsx';
    'rhox';
    'x0';
    'lambda';
    'xi'};

dispNames = {'$ M $';
    '$ \theta_{aux} $';
    %'$ \xi $';
    '$fixcost$';
    '$ \kappa $';
    '$ \zeta $';
    %'$ \psi $';
    '$ \epsilon_x $';
    '$ \rho_x $';
    '$ \bar{x} $';
    '$ \lambda $';
    '$ x_0 $'};

description = {'mass of potential entrants';
    'resale value of owned capital (aux)';
    %'share of rented capital';
    'addition fixed cost';
    'total capital per small firm';
    'marginal utility of leisure';
    %'exogenous firm exit rate';
    'standard deviation of log(x)';
    'autocorrelation of log(x)';
    'mean of log(x)';
    'initial debt distrib.';
    'productivity shifter of entrants'};

% Cell N*3 of characters for Table for exo parameters
% col1: field name; col2: name for Latex; col3: description
ExoNames = ...
    {'beta'   , '$\beta$',       'Subjective discount factor';
    'alpha'   , '$\alpha$',      'Capital Share';
    'delta_k' , '$\delta_k$',    'Capital depreciation rate 1';
    'delta_ks', '$\delta_{ks}$', 'Capital depreciation rate 2';
    'gamma'   , '$\gamma$',      'Span of control';
    'A'       , '$A$',           'TFP shifter';
    'psi'     , '$\psi$',        'Exogenous exit rate';
    'fixp'    , '$\xi$',         'Share of rented capital'};
%'theta_aux','$\theta$',      'Resale value of owned capital'


if ~isequal(numel(calibNames),numel(description))
    error("character arrays <calibNames> and <description> MUST have the same number of elements")
end

if ~isequal(numel(calibNames),numel(dispNames))
    error("character arrays <calibNames> and <dispNames> MUST have the same number of elements")
end

% guess must be a column vector
guess     = struct2vec(par,calibNames);
par.theta = par.theta_aux*(1-par.fixp);


%% Load data moments for the steady state
% Source: shared forlder, data, moments, data_moments.txt
% Remove 'autocorr_emp'
targetNames = ...
    {'avefirmsize';
    'empshare_small';
    'revshare_small';
    'exitrate';
    %'exitrate_age0';
    'avefirmsize_age0';
    'avefirmsize_age1';
    'jcr';
    'jdr';
    'rent_to_fixedcost';
    'fixedcost_to_rev';
    'fixedcost_to_payroll';
    'autocorr_emp';
    'debt_asset_all';
    'debt_asset_entrants';
    'ave_work';
    'hasNetDebt';
    %'hasNetDebt_age0';
    'debt_payroll_cond';
    'cash_payroll_cond';
    'k_payroll_ratio'};
% Moments not added to the above list:
% autocorr_logrev,fixedcost_to_payroll


targetNames_long = {'Average employment in small firms';
    'Small firm share of employment';
    'Small firm share of revenues';
    'Small firm exit rate';
    %'Small firm exit rate, age 0';
    'Average employment, age 0';
    'Average employment, age 1';
    'Job creation rate';
    'Job destruction rate';
    'Rent to fixed cost ratio';
    'Fixed expense to revenue ratio';
    'Fixed expense to payroll ratio';
    'Autocorr. employment';
    'Debt to asset ratio';
    'Debt to asset ratio of entrants';
    'Time spent in market work';
    'Share of firms with debt';
    %'Share of entrans with debt';
    'Debt to payroll cond';
    'Cash to payroll cond';
    'Capital payroll ratio'};

if ~isequal(numel(targetNames),numel(targetNames_long))
    error("character arrays <targetNames> and <targetNames_long> MUST have the same number of elements")
end

calibWeights = ones(numel(targetNames),1);
calibWeights(1) = 20; % avefirmsize

calibWeights(3) = 0; % revshare_small
calibWeights(4) = 20; % exitrate
%calibWeights(5) = 0; % exitrate_age0

calibWeights(5) = 8; % avefirmsize_age0
calibWeights(6) = 0; % avefirmsize_age1
calibWeights(7) = 1; % jcr

calibWeights(8) = 0; % jdr
calibWeights(9) = 0; % rent_to_fixedcost
calibWeights(11) = 0; %fixedcost_to_payroll
calibWeights(12) = 0.5; % autocorr_emp
calibWeights(13) = 1; % debt_asset_all
calibWeights(14) = 1; % debt_asset_entrants
calibWeights(15) = 1; % ave_work

calibWeights(16) = 1; % hasNetDebt
calibWeights(17) = 0; % debt_payroll_cond
calibWeights(18) = 0; % cash_payroll_cond
calibWeights(19) = 0.0; % k_payroll_ratio
filename = fullfile('..','data_moments','data_moments.txt');
data_mom = read_data_targets(filename); %data_mom is a structure



%% Comparative statics

if par.do_calib == 0
    disp('Compute steady-state at given parameter values')
    [obj_smm,sol,agg,b_grid,distribS,prices,model_mom,flag_ss,par] = ...
        fun_obj(guess,par,bounds,calibNames,data_mom,targetNames,...
        calibWeights,description,dispNames,targetNames_long);
    if flag_ss<0
        warning('Steady-state computation failed!')
    else
        if agg.K_corp<0
            warning("Capital in corporate sector is negative!")
        end
    end
    
else
    
    disp('Start comparative statics..')
    
    % Comparative statics wrt rhox
    % Compute moments for each value of the parameter. Save the moments in
    % matrix mom_mat with dim:(nvary,nmom)
    
    parName   = 'mass';
    parLimits = [0.01,4.45];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'theta_aux';
    parLimits = [0.5,0.99];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'fixcost';
    parLimits = [0.4,0.7];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'kappa';
    parLimits = [10,70];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'epsx';
    parLimits = [0.02,0.04];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'rhox';
    parLimits = [0.94,0.973];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'x0';
    parLimits = [3.5,4];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'lambda';
    parLimits = [0.01,0.05];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
    
    parName   = 'xi';
    parLimits = [0.3,0.9];
    [mom_all,mom_mat] = fun_perturb_param(parName,parLimits,targetNames,par);
  
end %END par.do_calib

% Write nvary*nmom matrix "mom_all" to the Excel file "calibration.xlsx"
% skipping the first row
% First row of the excel sheet must be 


