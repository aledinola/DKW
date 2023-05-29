%% Di Nola, Kaas, Wang (2022)
% FOLDER: 'alternative'
%{
----------------------- LEGEND --------------------------------------------

Rescue Policies for Small Businesses During the Covid-19 Recession
August 2022.
This main script calls <fun_steady_state> to compute the pre-pandemic
steady-state of the model. Then it calls <fun_transition> which computes
the transition of the economy after a one-period pandemic shock hits the
economy in t=1. Recall that t=1 in the code corresponds to t=0 in the
draft.
FIGURES:
    To make plots, run the following scripts: 
    "main_plots", "make_plots_compare","cum_impact_compare"
    This script loads the mat files saved in subfolder mat, draws figures
    and saves the figures in subfolder figures. The figures are then
    imported in the latex document in subfolder writeup.
TABLES:
    To make tables, run the following script: "main_tables". 
    This script loads saved mat files, generates tables in subfolder "tables". 
    The tables are read by the tex file in "writeup".

% 1 = impact, grant,  2 = no impact, grant, 3 = impact, no grant, 4 =  no impact, no grant

----------------------- END LEGEND ----------------------------------------
%}

clear
clc
close all
format long g
digits(100)
% Add path to numerical tools
addpath(genpath(fullfile('tools')));

disp("Alternative calibration of productivity process")

global obj_smm_best obj_tran_best %#ok<GVMIS> 


%% Set flags

par.do_calib = 0; % 0 = steady-state,
% 1 = steady-state+transition,
% 2 = calibration steady-state,
% 3 = calibration of transition shocks

par.grant_flag =0; % 0 = no grant; 1 = baseline grant (uniform); 2 = targeted grant by size
% 3 = uniform grant large; 4 = uniform grant small
% 5 = eta_g==eta_i (only impacted firms get grant)
% 6 = eta_g==eta_i (only impacted firms get grant), large grant
% 7 = eta_g==eta_i (only impacted firms get grant), small grant

par.grant_target = 0; % 0 = grant is untargeted (baseline); 
%                       1 = grant is targeted to impacted firms, some unimpacted firms also get grant;
%                       2 = slim targeted: only impacted firms get grant.

par.FigDir   = fullfile('figures','grant'); % folder to save figures
par.TabDir   = fullfile('tables'); % folder to save tables
par.InpDir   = fullfile('inputs'); % folder to read parameters
do_plots     = 0; %flag 0/1 to draw plots
par.do_plots_ss = 0; % flag 0/1 to draw plots for steady-state
do_save      = 1; %flag 0/1 to save plots as png and mat files
par.do_table = 1; %flag 0/1 to write tables on screen
par.do_tex   = 0; %flag 0/1 to write tables on tex (only if do_table=1)
par.verbose  = 1; %flag 0/1/2: 0 no display at all, 1=moderate, 2=disp everything
par.disp_mu  = 0;     % flag 0/1 to display iterations of mu
par.disp_tran = 1; %flag 0/1
do_debug     = 0;
est_algo     = 'fminsearch'; %Available options: 'simulan','simulannealbnd','fminsearch'
makeCompleteLatexDocument = 0; %flag 0/1 to generate a stand-alone tex doc
file_params  = 'estim_params.txt'; % Name of txt file where estimated parameters are saved
file_shocks  = 'estim_shocks.txt'; % Name of txt file where estimated shocks are saved
% Check if folders exist and if not create them
if ~(isfolder(par.FigDir))
    disp("Folder does not exist, creating it now..")
    mkdir(par.FigDir)
end
if ~(isfolder(par.TabDir))
    disp("Folder does not exist, creating it now..")
    mkdir(par.TabDir)
end

%% Set parameters, initial guess for estimation, bounds and exogenous grids
[par,guess,bounds,calibNames,dispNames,description,ExoNames] = set_parameters(par,file_params);

%% Load data moments for the steady state, set calibration weights
[targetNames,targetNames_long,calibWeights,data_mom] = set_targets_ss();

%% Set shocks
[par,bounds_shocks,data_mom_trans,calibWeightsTran] = set_shocks(par,file_shocks);

%% Set grant and rescue policies parameters Xp and eta for transition

[par] = set_grant(par);

%% Call fun_obj

if par.do_calib == 0
    disp("Run the model (steady-state) at given parameter values")
    tic
    [obj_smm,sol,agg,b_grid,distribS,prices,model_mom,flag_ss,par] = fun_obj(guess,par,bounds,calibNames,data_mom,targetNames,calibWeights,description,dispNames,targetNames_long);
    toc
    % Save workspace in subfolder \mat
    if do_save==1
        disp('Saving steady-state results as mat file..')
        save(fullfile('mat','ss.mat'))
        disp('Done!')
    end
    
elseif par.do_calib == 1
    disp("Steady-state and transition")
    
    % Solve for the steady-state
    fprintf('Start steady-state computation... \n')
    fprintf(' \n')
    tic
    [sol,agg,b_grid,distribS,prices,model_mom,flag_ss,par] = fun_steady_state(par);
    if flag_ss<0
        warning("steady-state not found!")
    else
        fprintf('Steady-state computed! \n')
    end
    toc
    
    fprintf(' \n')
    fprintf('Start transition... \n')
    fprintf('Transition periods T       = %d \n',par.T)
    fprintf('Frac. impacted firms eta_i = %f \n',par.eta_i)
    fprintf('grant_flag                 = %d \n',par.grant_flag)
    fprintf('grant_target               = %d \n',par.grant_target)
    %fprintf('Xp                         = %f \n',par.Xp)
    fprintf(' \n')

    tic
    [agg_tran,path,conv_flag,pol_tran,distrib_tran] = fun_transition(par,sol,agg,distribS,prices,b_grid);
    disp('Time to compute transition:')
    toc
    if conv_flag<0
        warning("Be careful!")
        fprintf('Transition did not converge \n' )
    end
    
    % Compute moments for transition
    %if conv_flag>=0
    [model_mom_trans,irf] = fun_targets_tran(data_mom_trans,agg_tran,path,agg,prices,calibWeightsTran);
    %end
    
    % Print shock parameters to screen
    fprintf("eta_i:         %8.6f  \n",par.eta_i);
    fprintf("v_small:       %8.6f  \n",par.v_small);
    fprintf("v_corp:        %8.6f  \n",par.v_corp);
    fprintf("util_shift:    %8.6f  \n",par.util_shift);
    fprintf("lsupply_shift: %8.6f  \n",par.lsupply_shift);
    fprintf("mass_shift: %8.6f  \n",par.mass_shift);
    fprintf("lambda_shift: %8.6f  \n",par.lambda_shift);
    fprintf("rho_shock:     %8.6f  \n",par.rho_shock);
    fprintf(" ")
    % Save transition results
    results.agg_tran = agg_tran;
    results.path     = path;
    
    % Print calibration transition to txt file
    txt_export_tran(par,model_mom_trans,agg_tran,par.TabDir);


    % Save workspace in subfolder \mat
    if do_save==1
        pol_tran_small.pol_exit   = pol_tran.pol_exit;
        pol_tran_small.b_grid     = pol_tran.b_grid; 
        pol_tran_small.V1         = pol_tran.V1; % needed for zombie firms
        pol_tran_small.pol_kp_unc = pol_tran.pol_kp_unc; % needed for grant induced overinvestment

        pol_tran = pol_tran_small;
        clearvars -except par agg agg_tran b_grid data_mom data_mom_trans ...
            prices path model_mom_trans data_mom_trans irf pol_tran ...
            distrib_tran sol do_plots do_save

        if par.grant_flag==1 && par.grant_target == 0
            save(fullfile('mat','grant_baseline.mat'),'-v7.3')
            SaveDir = fullfile('figures','grant_baseline');
        elseif par.grant_target == 1
            save(fullfile('mat','grant_targeted.mat'),'-v7.3')
        elseif par.grant_target == 2 && par.grant_flag == 5
            save(fullfile('mat','grant_targslim.mat'),'-v7.3')
            SaveDir = fullfile('figures','grant_targslim');
        elseif par.grant_target == 2 && par.grant_flag == 6
            save(fullfile('mat','grant_targslim_large.mat'),'-v7.3')
            SaveDir = fullfile('figures','grant_targslim_large');
        elseif par.grant_target == 2 && par.grant_flag == 7
            save(fullfile('mat','grant_targslim_small.mat'),'-v7.3')
            SaveDir = fullfile('figures','grant_targslim_small');
        elseif par.grant_flag==5 && par.grant_target == 0
            save(fullfile('mat','grant_eta.mat'),'-v7.3')
        elseif par.grant_flag==0
            save(fullfile('mat','nogrant.mat'),'-v7.3')
            SaveDir = fullfile('figures','nogrant');
        elseif par.grant_flag==3 && par.grant_target == 0
            save(fullfile('mat','grant_large.mat'),'-v7.3')
        elseif par.grant_flag==4 && par.grant_target == 0
            save(fullfile('mat','grant_small.mat'),'-v7.3')
        end
        
        % Plot IRFs of transition and save them
        plot_irf(irf,SaveDir,do_save,'-dpng')

    end
    
elseif par.do_calib == 2
    disp("Start calibration of the steady-state..")
    
    % Delete txt file results_sofar if present
    if isfile(fullfile('inputs','results_sofar.txt'))
        delete(fullfile('inputs','results_sofar.txt'));
    end
    
    obj_smm_best = realmax;
    
    % Define the function to be minimized
    f_obj = @(x) fun_obj(x,par,bounds,calibNames,data_mom,targetNames,calibWeights,description,dispNames,targetNames_long);
    %Convert bounds from struct to vector
    bounds_vec = bounds2vec(bounds,calibNames);
    lbounds_vec = bounds_vec(:,1); %LOWER bounds
    ubounds_vec = bounds_vec(:,2); %UPPER bounds
    
    switch est_algo
        
        case 'simulan'
            disp('Simulated annealing')
            %%%%%%%%%%%%%% SIMULATED ANNEALING
            %%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            maxim = 0;       % 1 if maximize, 0 otherwise
            rt    = 0.85;    % temperature decreasing ratio
            eps   = 1e-4;    % error tolerance for termination
            ns    = 20;      % number of cycles
            nt    = 100; %max(100,5*length(guess)); % number of iterations before temperature reduction
            neps  = 4;       % number of final functions values used to decide upon termination
            maxevl= 100000;  % maximum number of function evaluations
            iprint= 1;%2;    % control inside printing
            t     = 0.5;     % initial temperature
            option=[maxim rt eps ns nt neps maxevl iprint t]';
            c     = 2*ones(length(guess),1);
            vm    = 0.1*ones(length(guess),1);
            
            [x,fopt,nacc,nfcnev,nobds,ier,t,vm]= SIMULANS(f_obj,guess,option,lbounds_vec,ubounds_vec,c,vm);
            disp('****   results after sa   ****');
            disp('solution');
            disp(x(:)');
            disp('final step length');
            disp(vm(:)');
            fprintf('optimal function value              : %g\n',fopt);
            fprintf('number of function evaluations      : %g\n',nfcnev);
            fprintf('number of accepted evaluations      : %g\n',nacc);
            fprintf('number of out-of-bounds evaluations : %g\n',nobds);
            fprintf('final temperature                   : %g\n',t);
            fprintf('error                               : %g\n',ier);
            
        case 'simulannealbnd'
            % Similar to simulan, but uses a Matlab-provided routine. See
            % Matlab documentation on simulannealbnd.
            disp('Simulated annealing: simulannealbnd')
            init_temp = 5;
            options = optimoptions(@simulannealbnd,'Display','iter','InitialTemperature',init_temp);
            [x,~,exitflag] = simulannealbnd(f_obj,guess,lbounds_vec,ubounds_vec,options);
            if exitflag<=0
                warning('simulannealbnd did not converge')
            end
            
        case 'fminsearch'
            disp('Nelder-Mead')
            options = optimset('MaxIter',5000,'Display','iter');
            %[x,~,exitflag] = fminsearch(f_obj,guess,options);
            [x,~,exitflag] = fminsearchcon(f_obj,guess,lbounds_vec,ubounds_vec,[],[],[],options);
            if exitflag<=0
                warning('fminsearch did not converge')
            end
            
        otherwise
            error("Selected algo for estimation is not available")
    end
    
    % Save estimated parameters in file estim_params.txt
    FID = fopen(fullfile(par.InpDir,'estim_params.txt'),'w');
    for i = 1:length(x)
        fprintf(FID,'%s  %.10f \n',calibNames{i},x(i));
    end
    fclose(FID);
    
elseif par.do_calib == 3
    disp("START CALIBRATION OF TRANSITION")
    % 6 parameters (shock values)
    % Many targets (delta GDP,C,I,output_small)
    
    % Delete txt file results_sofar if present
    if isfile(fullfile('inputs','results_tran_sofar.txt'))
        delete(fullfile('inputs','results_tran_sofar.txt'));
    end
    
    obj_tran_best = realmax;
    
    options = optimset('MaxIter',500,'Display','iter');
    f = @(x_in) fun_calib_transition(x_in,par,data_mom_trans,calibWeightsTran,bounds_shocks);
    X0 = [par.eta_i, v_corp, util_shift, lsupply_shift, rho_shock]';
    [x,~,exitflag] = fminsearch(f,X0,options);
    
    if exitflag<=0
        warning('fminsearch did not converge')
    end
    disp("CALIBRATION OF TRANSITION IS OVER")
    
    % Save estimated shocks in file estim_shocks.txt
    FID = fopen(fullfile(par.InpDir,'estim_shocks.txt'),'w');
    fprintf(FID,"eta_i:         %8.6f  \n",x(1));
    fprintf(FID,"v_corp:        %8.6f  \n",x(2));
    fprintf(FID,"util_shift:    %8.6f  \n",x(3));
    fprintf(FID,"lsupply_shift: %8.6f  \n",x(4));
    fprintf(FID,"rho_shock:     %8.6f  \n",x(5));
    fclose(FID);    
else
    error("Selected value for do_calib is not available")
    
end
