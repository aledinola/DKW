function [par,guess,bounds,calibNames,dispNames,description,ExoNames] = set_parameters(par,file_params)
% SET NUMERICAL AND ECONOMIC PARAMETERS OF THE MODEL HERE
% INPUTS:
%   "par"         Structure with model parameters
%   "file_params" Character variable with file name of input parameters
% OUTPUTS:
%   "par"         Updated structure
%   "guess"       Column vector with initial guess for calibration
%   "bounds"      Structure where each field is a 1*2 numeric vector
%   "calibNames"  Cell array of characters
%   "dispNames"   Cell array of characters
%   "description" Cell array of characters
%   "ExoNames"    Cell array of characters

%% Inputs check
if ~isstruct(par)
    error('input argument "par" must be a structure')
end
if ~ischar(file_params)
    error('input argument "file_params" must be a character')
end

%% Numerical parameters

par.T            = 180;  % length transition
par.Tmax         = 8;    % after Tmax, shocks are all zero
par.nx           = 60;   % grid size for productivity
par.nb           = 80;   % grid size for debt
par.nk           = 100;  % grid size for capital
par.ns           = 2;    % grant vs no grant
par.ni           = 2;    % impacted vs unimpacted
par.nn           = par.ns*par.ni;
par.x_process    = 1;   % 1 = AR1; 2 = bounded pareto distribution with persistence
par.k_distrib    = 2;   % 1 = uniform distribution; 2 = pareto distribution; 
par.k_lb         = 0.1; % lower bound of k_grid
par.k_ub         = 200;   % upper bound of k_grid
par.k_min        = par.k_lb; % minimum of the uniform distribution of capital for entrants
par.tol_bhat     = 1e-9; % Tolerance for fixed point B_hat(k,x)
par.tol_vfi      = 1e-9; % tolerance for VFI
par.max_iter     = 4000;  % max iter for VFI
par.tol_vfi_u    = 1e-9;  % tolerance for VFI of the unconstrained firms
par.do_howard    = 1;     % flag 0/1 Howard acceleration
par.n_howard     = 50;
par.tol_dist     = 1e-6; % tolerance for distribution
par.maxiter_dist = 10000; % max no iter for distribution
par.max_iter_tr  = 100;   % max no iter for transition
par.tol_tran     = 0.0005; % tol for transition
par.dampening    = 1;   % dampen update in bisection
par.seed         = 67354; % seed for random number generator, used in fun_simulate
par.write_calib_append = 0;
par.do_write_calib = 0;
par.N_sim = 80000; 
par.T_sim = 68; % 17 years
par.emp_min = 0; % threshold for microbusinesses. Microbusinesses are excluded from exit rate calculation.

%% Set exogenous parameters (not part of internal calibration)

par.beta      = 0.989;  % discount factor
par.sigma     = 2.0;    % CRRA parameter
par.alpha     = 0.3;    % Cobb-Douglas exponent on capital
par.delta_k   = 0.015;  % depreciation rate capital
par.gamma1    = 0.3182; % share of capital prod function small firms
par.gamma2    = 0.88;   % span of control parameter in f(l) small firms
par.A         = 0.25;   % production function shifter for the corporate sector
%par.psi       = 0.004; % exogenous exit rate
par.lambda0   = 1;      % Tightness of collateral constraint (i.e. lam=lam0*theta*(1-delta))
par.cost_e    = 0;
par.bk0_vec   = [-0.09375,0.125,0.8684]';  % initial debt to capital ratio (p25,p50,p75 of entrants, KFS)
par.bk0_prob  = [0.25,0.5,0.25]';  % probability for bk0
% productivity process parameters: AR1 process
par.x0        = 1;    % Ln(x0) mean of AR(1)
par.xi        = 1;    % productivity gap of potential entrants relative to incumbents
par.epsx      = 0.12;   % standard deviation of innovation in AR(1) for x'|x
par.rhox      = 0.95;   % persistence in AR(1) for x'|x

% Productivity process parameters: bounded pareto distribution for x
par.x_lb      =  0.5;
par.x_ub      =  4.0;
par.x_shape   =  0.1;
par.x_rho     =  0.9306;

% Initial capital distribution parameters: AR1 process
par.k_max = 41;
% Initial capital distribution parameters: Pareto distribution
par.k_alpha = 0.3240858974;

%% Initial conditions for internal parameters
% Read parameters from 'estim_params.txt'

fprintf('Reading parameters from file "%s" \n',file_params)
FID = fopen(fullfile(par.InpDir,file_params));
C   = textscan(FID,'%s %f');
fclose(FID);

names  = C{1};
values = C{2};

for i=1:numel(names)
    par.(names{i}) = values(i);
end

%% Set bounds for internal parameters
% [First number is lower bound, second number is the upper bound]

bounds.mass      = [0.0001, 10000];  % mass of potential entrants
bounds.fixcost1  = [0.0001, 1];      % additional fixed operation cost
bounds.fixcost2  = [0.000, 1];       % additional fixed operation cost
bounds.theta     = [0.1,1];          % Resale value of capital
%bounds.lambda0     = [1,1];         % lambda=lambda0*theta*(1-delta)
bounds.psi     = [0,0.01];           % exogenous exit rate
%bounds.k_min     = [1e-6,50];       % Pareto distrib for capital, lower bound
bounds.k_alpha   = [0.1,2];          % Pareto distrib for capital, shape
%bounds.k_max     = [1e-6,200];      % Uniform distribution for capital, upper bound
bounds.zeta      = [0.0001,10000.0]; % utility of leisure
%bounds.b0        = [-100, 1000];    % Initial debt level for entrants
bounds.x0        = [0.1, 10.0];      % Mean of log(x), productivity shock 
bounds.epsx      = [0.01,1];         % Dispersion of log(x)
bounds.rhox      = [0.5, 0.999];     % Persistence of log(x)     
% Set parameter names IN THE SAME ORDER AS THEY APPEAR in bounds and guess
% calibNames must be column cell array of characters

calibNames = {'mass';
    'fixcost1';
    'fixcost2';
    'theta';
    %'lambda0';
    'psi';
    %'k_min';
    'k_alpha';
    %'k_max';
    'x0';
    'epsx';
    'rhox';
    %'xi';
    'zeta'};
    %'b0'};

% For latex tables
dispNames = {'$ M $';
    '$ fixcost1  $';
    '$ fixcost2  $';
    '$ \theta $';
    '$ \psi $';
    '$ \alpha_{\kappa}  $';
    '$ \bar{x}    $';
    '$ \varepsilon_x    $';
    '$ \rho_x    $';
    '$ \zeta    $'};

description = {'mass of potential entrants';
    'Intercept fixed cost';
    'Slope fixed cost';
    'Resale value of capital';
    'Exogenous exit rate';
    'Shape Pareto capital';
    'Ln(x0) mean of AR(1)';
    'Dispersion of $x$';
    'Persistence of $x$';
    'marginal utility of leisure'};

% Cell N*3 of characters for Table for exo parameters
% col1: field name; col2: name for Latex; col3: description
ExoNames = ...
    {'beta'   , '$\beta$',       'Subjective discount factor';
    'sigma'   , '$\sigma$',      'CRRA coefficient';
    'alpha'   , '$\alpha$',      'Capital Share corporate sector';
    'delta_k' , '$\delta_k$',    'Capital depreciation rate';
    'lambda0' , '$\lambda_0$',    'Collateral constraint parameter';
    'gamma1'  , '$\gamma_1$',    'Capital Share small firms';
    'gamma2'  , '$\gamma_2$',    'Span of control';
    'A'       , '$A$',           'TFP shifter'};


if ~isequal(numel(calibNames),numel(description))
    error("character arrays <calibNames> and <description> MUST have the same number of elements")
end

if ~isequal(numel(calibNames),numel(dispNames))
    error("character arrays <calibNames> and <dispNames> MUST have the same number of elements")
end

% guess must be a column vector
guess     = struct2vec(par,calibNames);

%% Generate initial grid for idios productivity x

if par.x_process == 1
    %-- AR1, this is the defeault case
    par.mean_x             = (1-par.rhox)*log(par.x0);
    par.cover              = 3.5; % standard value, see Tauchen (1986)
    [par.pi_x, log_x_grid] = markovapprox(par.rhox,par.epsx,par.mean_x,par.cover,par.nx);
    par.x_grid             = exp(log_x_grid);
elseif par.x_process == 2
    %-- Bounded Pareto with persistence
    par.x_grid = linspace(par.x_lb,par.x_ub,par.nx);
    [par.pi_x, par.x_prob] = paretojo(par.nx,par.x_grid,par.x_shape,par.x_rho);

else
    error("set_parameters: x_process is invalid!")
end

%% Generate grid for capital

par.k_grid = linspace(par.k_lb,par.k_ub,par.nk)';


end %end function "set_parameters"

