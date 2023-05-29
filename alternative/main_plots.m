%% This script generates figures (and some tables/other results) for the paper
% The script loads the following mat files:
% ss.mat, grant_baseline.mat, nogrant.mat, grant_targslim.mat
% which must be saved in subfolder "mat"
% - Plot steady-state distributions, policy functions (exit, entry, investment)
% - Plot calibration of transition, and plot the exit thresholds in t=1 for the baseline grant
% - Plot change in exit rate due to grant in t=1
% - Decomposing output change in the pandemic
% - Plot reallocation among small firms (across bins of x)
% - Four types of capital adjustment
% - Compute fraction of saved firms and zombie firms
% - Plot four types of capital adjustment
% - Plots ave b and ave b/k in transition, grant vs no grant
% - Welfare analysis

clear;clc;close all
% Add path to numerical tools
addpath(genpath(fullfile('tools')));
% %% Set up some useful paths
% ResultsDir = fullfile('mat'); %folder where .mat files are stored
FormatFig = '-depsc';  % Specify '-dpng' or '-depsc'

matNames = {'ss','nogrant','grant_baseline','grant_targslim'};
for ii = 1:numel(matNames)
    if ~isfile(fullfile('mat',[matNames{ii},'.mat']))
        warning('MAT file ''%s'' is missing \n',matNames{ii})
        pause
    end
end

%% Plot Steady-state distributions, policy functions (exit, entry, investment)
disp('Plot Steady-state distributions, policy functions')
% Load results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'b_grid','sol','par','prices','distribS','model_mom','data_mom')
else
    error('File ss.mat does not exist!')
end

% Specify here where you want to save the figures
SaveDir = fullfile('figures','ss');

FS  = 16; % font size for plots

% Call function to make plots and to save them
plot_ss(b_grid,sol,par,prices,distribS,model_mom,data_mom,SaveDir,1,FormatFig,FS);

% We plot policy function of investment
plot_ss_policy(b_grid,sol,distribS,par,SaveDir,FormatFig,1,FS)

%% Plot calibration of transition, and plot the exit thresholds in t=1 for the baseline grant
clear; close all
disp('Plot calibration of transition, and plot the exit thresholds')
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'model_mom_trans','par','b_grid','sol','pol_tran');
else
    error('File grant_baseline.mat does not exist')
end

FormatFig = '-depsc';  % Specify '-dpng' or '-depsc'
FS        = 16; % font size for plots
do_save   = 1;
SaveDir   = fullfile('figures','grant_vs_nogrant');
% Calibration of transition: moments and shocks
plot_calib_tran(model_mom_trans,par,FormatFig,FS,do_save,SaveDir);
% Exit thresholds: grant receiving and no-grant receiving
[exit_rule_ss,exit_rule_grant_baseline] = plot_exit_tran(par,b_grid,sol,pol_tran,SaveDir,do_save,FormatFig);

%% Financial saving and debt along the transition
disp('Financial saving and debt along the transition')
clear;close all

% Steady-state
load(fullfile('mat','ss.mat'),'par','b_grid','distribS');
nb = par.nb;
k_grid = par.k_grid;
mu_active = distribS.mu_active;
ave_b_ss = sum(b_grid.*mu_active,'all')/sum(mu_active,'all');
bk_mat = b_grid./repmat(k_grid,1,nb); %(nk,nb)
ave_bk_ss = sum(bk_mat.*mu_active,'all')/sum(mu_active,'all');

% No grant economy
load(fullfile('mat','nogrant.mat'),'par','pol_tran','distrib_tran');
b_grid_nogrant       = pol_tran.b_grid;
distrib_tran_nogrant = distrib_tran;
weights_nogrant      = par.weights;
% outputs have dim (T+1,1)
[ave_b_nogrant,ave_bk_nogrant] = plot_liquidity(b_grid_nogrant,distrib_tran_nogrant,weights_nogrant,par);

% Baseline grant economy
load(fullfile('mat','grant_baseline.mat'),'par','pol_tran','distrib_tran');
b_grid_grant_baseline       = pol_tran.b_grid;
distrib_tran_grant_baseline = distrib_tran;
weights_grant_baseline      = par.weights;
% outputs have dim (T+1,1)
[ave_b_grant_baseline,ave_bk_grant_baseline] = plot_liquidity(b_grid_grant_baseline,distrib_tran_grant_baseline,weights_grant_baseline,par);


% Percentage deviation from steady-state
irf_ave_b_grant_baseline = 100*(ave_b_grant_baseline-ave_b_ss)./abs(ave_b_ss);
irf_ave_b_nogrant = 100*(ave_b_nogrant-ave_b_ss)./abs(ave_b_ss);

irf_ave_bk_grant_baseline = 100*(ave_bk_grant_baseline-ave_bk_ss)./abs(ave_bk_ss);
irf_ave_bk_nogrant = 100*(ave_bk_nogrant-ave_bk_ss)./abs(ave_bk_ss);

% Plot options
FormatFig = '-dpng';  % Specify '-dpng' or '-depsc'
FS = 16; % font size for plots
LW = 3; % line width
T_last = min(180,par.T);
do_save = 1;
SaveDir = fullfile('figures','grant_vs_nogrant');

% Do plots

% Average b, baseline grant vs no grant
figure
plot(0:T_last-1,irf_ave_b_grant_baseline(1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,irf_ave_b_nogrant(1:T_last),'--','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Baseline grant','Laissez-faire','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'irf_ave_b'),FormatFig); end

% Average b/k, baseline grant vs no grant
figure
plot(0:T_last-1,irf_ave_bk_grant_baseline(1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,irf_ave_bk_nogrant(1:T_last),'--','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Baseline grant','Laissez-faire','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'irf_ave_bk'),FormatFig); end


%% Plot change in exit rate due to grant in t=1
clear;close all
disp('Plot change in exit rate due to grant in t=1')
% Load results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'sol','prices','distribS')
else
    error('File ss.mat does not exist!')
end
% No grant or laissez-faire economy
if isfile(fullfile('mat','nogrant.mat'))
    load(fullfile('mat','nogrant.mat'),'par','pol_tran');
else
    error('File nogrant.mat does not exist')
end
pol_tran_nogrant     = pol_tran;
weights_nogrant      = par.weights;

% Baseline grant economy
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'par','pol_tran');
else
    error('File grant_baseline.mat does not exist')
end
pol_tran_grant_baseline     = pol_tran;
weights_grant_baseline      = par.weights;

% Compute statistics by employment l
%perct = linspace(1/4,3/4,3)';
%nbins = length(perct)+1;

nbins = 4;
% TODO: make the bins as follows
% [0-10, 10-20, 20-100, 100+]
[exitSS,exit_nogrant,exit_grant_baseline,binlids] = plot_micro_l(sol,distribS,...
    prices,pol_tran_nogrant,weights_nogrant,pol_tran_grant_baseline,...
    weights_grant_baseline,par);

% Change in exit rate from no grant economy to grant economy
change_exit_rate  = exit_grant_baseline.exit_rate - exit_nogrant.exit_rate;
change_exit       = exit_grant_baseline.exits - exit_nogrant.exits;
share_change_exit = change_exit / sum(change_exit);

% Options for plots
FS        = 16; % font size
FORMAT    = '-depsc'; % Specify '-dpng' or '-depsc'
SaveDir   = fullfile('figures','grant_vs_nogrant');
do_save   = 1;
bin_names = {'0-9';'10-19';'20-99';'100+'};

% Employment share in ss vs share of change in exit due to grant, by decile of x
figure
yyaxis left
w1 = 0.8;
bar(1:nbins,exitSS.emp_share,w1)
%text(1:nbins, aveemp_bin_SS/10, num2str(aveemp_bin_SS,'%6.2f'),'vert','bottom','horiz','center')
hold on
yyaxis right
w2 = 0.6;
bar(1:nbins,change_exit_rate,w2)
set(gca,'FontSize',12);
set(gca,'xticklabel',bin_names)
xlabel('Bins by firm size (employment)','FontSize',FS)
%title('Effect of grant on exit rate','FontSize',14)
legend('Employment share in SS','Change in exit rate','FontSize',FS,'location','southoutside')
hold off
if do_save==1; print(fullfile(SaveDir,'emp_exit_rate_lbins'),FORMAT); end

%% Decomposing output change in the pandemic
disp('Decomposing output change in the pandemic')
clear; close all
% Three factors: (1) TFP, (2) Change in employment policy, and (3) Change in entry, exit, and investment
% Load results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'prices','distribS')
else
    error('File ss.mat does not exist!')
end
% No grant or laissez-faire economy
if isfile(fullfile('mat','nogrant.mat'))
    load(fullfile('mat','nogrant.mat'),'par','distrib_tran','path');
else
    error('File nogrant.mat does not exist')
end
%pol_tran_nogrant     = pol_tran;
distrib_tran_nogrant = distrib_tran;
path_nogrant         = path;
weights_nogrant      = par.weights;

% Baseline grant economy
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'par','distrib_tran','path');
else
    error('File grant_baseline.mat does not exist')
end
%pol_tran_grant_baseline     = pol_tran;
distrib_tran_grant_baseline = distrib_tran;
path_grant_baseline         = path;
weights_grant_baseline      = par.weights;

% Slim targeted grant
if isfile(fullfile('mat','grant_targslim.mat'))
    load(fullfile('mat','grant_targslim.mat'),'par','distrib_tran','path');
else
    error('File grant_targslim.mat does not exist')
end
%pol_tran_grant_targslim     = pol_tran;
distrib_tran_grant_targslim = distrib_tran;
path_grant_targslim         = path;
weights_grant_targslim      = par.weights;

T_last = 8;
Delta_grant_baseline = fun_decomposition(T_last,distribS,prices,distrib_tran_grant_baseline,path_grant_baseline,weights_grant_baseline,par);
Delta_nogrant = fun_decomposition(T_last,distribS,prices,distrib_tran_nogrant,path_nogrant,weights_nogrant,par);
%Delta_grant_targeted = fun_decomposition(T_last,distribS,prices,distrib_tran_grant_targeted,path_grant_targeted,weights_grant_targeted,par);
Delta_grant_targslim = fun_decomposition(T_last,distribS,prices,distrib_tran_grant_targslim,path_grant_targslim,weights_grant_targslim,par);

% baseline
frac_deltaExit_grant_baseline = sum(Delta_grant_baseline.deltaExit,'all')/sum(Delta_grant_baseline.deltaY,'all');
frac_deltaTFP_grant_baseline = sum(Delta_grant_baseline.deltaTFP,'all')/sum(Delta_grant_baseline.deltaY,'all');
frac_deltaL_grant_baseline = sum(Delta_grant_baseline.deltaL,'all')/sum(Delta_grant_baseline.deltaY,'all');
% no grant
frac_deltaExit_nogrant = sum(Delta_nogrant.deltaExit,'all')/sum(Delta_nogrant.deltaY,'all');
frac_deltaTFP_nogrant = sum(Delta_nogrant.deltaTFP,'all')/sum(Delta_nogrant.deltaY,'all');
frac_deltaL_nogrant = sum(Delta_nogrant.deltaL,'all')/sum(Delta_nogrant.deltaY,'all');
% targeted grant
%frac_deltaExit_grant_targeted = sum(Delta_grant_targeted.deltaExit,'all')/sum(Delta_grant_targeted.deltaY,'all');
%frac_deltaTFP_grant_targeted = sum(Delta_grant_targeted.deltaTFP,'all')/sum(Delta_grant_targeted.deltaY,'all');
%frac_deltaL_grant_targeted = sum(Delta_grant_targeted.deltaL,'all')/sum(Delta_grant_targeted.deltaY,'all');
% slim targeted grant
frac_deltaExit_grant_targslim = sum(Delta_grant_targslim.deltaExit,'all')/sum(Delta_grant_targslim.deltaY,'all');
frac_deltaTFP_grant_targslim = sum(Delta_grant_targslim.deltaTFP,'all')/sum(Delta_grant_targslim.deltaY,'all');
frac_deltaL_grant_targslim = sum(Delta_grant_targslim.deltaL,'all')/sum(Delta_grant_targslim.deltaY,'all');


% Write to TXT file
FID = fopen(fullfile('tables','decomp_output_small.txt'),'w');
fprintf(FID,'--------------------------\n');
fprintf(FID,'Decomposition of output change from steady state \n');
fprintf(FID,'Number of periods: %d \n',T_last);
fprintf(FID,'Baseline Grant \n');
fprintf(FID,'Fraction of output change due to investment change or firm exit: %f \n',frac_deltaExit_grant_baseline);
fprintf(FID,'Fraction of output change due to TFP: %f \n',frac_deltaTFP_grant_baseline);
fprintf(FID,'Fraction of output change due to labor demand change: %f \n',frac_deltaL_grant_baseline);
fprintf(FID,'--------------------------\n');
fprintf(FID,'No Grant \n');
fprintf(FID,'Fraction of output change due to investment change or firm exit: %f \n',frac_deltaExit_nogrant);
fprintf(FID,'Fraction of output change due to TFP: %f \n',frac_deltaTFP_nogrant);
fprintf(FID,'Fraction of output change due to labor demand change: %f \n',frac_deltaL_nogrant);
fprintf(FID,'--------------------------\n');
fprintf(FID,'Slim targeted grant \n');
fprintf(FID,'Fraction of output change due to investment change or firm exit: %f \n',frac_deltaExit_grant_targslim);
fprintf(FID,'Fraction of output change due to TFP: %f \n',frac_deltaTFP_grant_targslim);
fprintf(FID,'Fraction of output change due to labor demand change: %f \n',frac_deltaL_grant_targslim);
fprintf(FID,'--------------------------\n');
fclose(FID);

% Write to TEX file
FID = fopen(fullfile('tables','decomp_output_small.tex'),'w');
fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
%fprintf(FID,'Decomposition of output change from steady state & & \n');
%fprintf(FID,'Number of periods: %d \n',T_last);
fprintf(FID,'Baseline Grant &  \\\\  \n');
fprintf(FID,'$ \\Delta_{K,Exit}/\\Delta_{Y} $ & %8.4f \\\\ \n',frac_deltaExit_grant_baseline);
fprintf(FID,'$\\Delta_{TFP}/\\Delta_{Y}  $   & %8.4f \\\\ \n',frac_deltaTFP_grant_baseline);
fprintf(FID,'$ \\Delta_{L}/\\Delta_{Y}    $   & %8.4f \\\\ \n',frac_deltaL_grant_baseline);
fprintf(FID, '\\hline \n');
fprintf(FID,'No Grant       &  \\\\ \n');
fprintf(FID,'$\\Delta_{K,Exit}/\\Delta_{Y} $ & %8.4f \\\\ \n',frac_deltaExit_nogrant);
fprintf(FID,'$\\Delta_{TFP}/\\Delta_{Y}   $  & %8.4f \\\\ \n',frac_deltaTFP_nogrant);
fprintf(FID,'$\\Delta_{L}/\\Delta_{Y}     $  & %8.4f \\\\ \n',frac_deltaL_nogrant);
fprintf(FID, '\\hline \n');
fprintf(FID,'Targeted grant &  \\\\ \n');
fprintf(FID,'$\\Delta_{K,Exit}/\\Delta_{Y}$  & %8.4f \\\\ \n',frac_deltaExit_grant_targslim);
fprintf(FID,'$\\Delta_{TFP}/\\Delta_{Y}   $  & %8.4f \\\\ \n',frac_deltaTFP_grant_targslim);
fprintf(FID,'$\\Delta_{L}/\\Delta_{Y}     $  & %8.4f \\\\ \n',frac_deltaL_grant_targslim);
fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);

%% Plots for reallocation among small firms (across bins of x)
%  small firm mass, employment, capital, and output by x in the short- and
%  long-run.
disp('Plots for reallocation among small firms (across bins of x)')
clear; close all
%
% Load results
load(fullfile('mat','ss.mat'),'distribS','prices')
% Baseline grant economy
load(fullfile('mat','grant_baseline.mat'),'distrib_tran','path','par')
distrib_tran_grant_baseline = distrib_tran;
path_grant_baseline         = path;
weights_grant_baseline      = par.weights;
% Laissez faire economy
load(fullfile('mat','nogrant.mat'),'distrib_tran','path','par')
distrib_tran_nogrant = distrib_tran;
path_nogrant         = path;
weights_nogrant      = par.weights;

% Compute small firm mass, employment, capital, and output by deciles of x and t
nbins = 2;
result_grant_baseline = plot_micro_t(nbins,distribS,prices,distrib_tran_grant_baseline,path_grant_baseline,weights_grant_baseline,par);
result_nogrant = plot_micro_t(nbins,distribS,prices,distrib_tran_nogrant,path_nogrant,weights_nogrant,par);

% Percentage deviation from steady-state
mass_dev_grant_baseline = 100*(result_grant_baseline.mass_xt-result_grant_baseline.mass_ss)./result_grant_baseline.mass_ss;
emp_dev_grant_baseline  = 100*(result_grant_baseline.emp_xt-result_grant_baseline.emp_ss)./result_grant_baseline.emp_ss;
K_dev_grant_baseline    = 100*(result_grant_baseline.K_xt-result_grant_baseline.K_ss)./result_grant_baseline.K_ss;

mass_dev_nogrant = 100*(result_nogrant.mass_xt-result_nogrant.mass_ss)./result_nogrant.mass_ss;
emp_dev_nogrant  = 100*(result_nogrant.emp_xt-result_nogrant.emp_ss)./result_nogrant.emp_ss;
K_dev_nogrant    = 100*(result_nogrant.K_xt-result_nogrant.K_ss)./result_nogrant.K_ss;

% Plot options
FormatFig = '-depsc'; % Specify '-dpng' or '-depsc'
FS      = 20;         % Font size
LW      = 4;          % Line width
T_last  = 80;
do_save = 1;
SaveDir = fullfile('figures','grant_vs_nogrant');

% Do plots

% mass, baseline grant
figure
plot(0:T_last-1,mass_dev_grant_baseline(1,1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,mass_dev_grant_baseline(2,1:T_last),'--','LineWidth',LW)
%hold on
%plot(0:T_last-1,mass_dev_grant_baseline(4,1:T_last),':','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Quartile 1','Quartile 2','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'mass_dev_grant_baseline'),FormatFig); end

% mass, no grant
figure
plot(0:T_last-1,mass_dev_nogrant(1,1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,mass_dev_nogrant(2,1:T_last),'--','LineWidth',LW)
%hold on
%plot(0:T_last-1,mass_dev_nogrant(4,1:T_last),':','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Quartile 1','Quartile 2','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'mass_dev_nogrant'),FormatFig); end

% employment, baseline grant
figure
plot(0:T_last-1,emp_dev_grant_baseline(1,1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,emp_dev_grant_baseline(2,1:T_last),'--','LineWidth',LW)
%hold on
%plot(0:T_last-1,emp_dev_grant_baseline(4,1:T_last),':','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Quartile 1','Quartile 2','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'emp_dev_grant_baseline'),FormatFig); end

% employment, no grant
figure
plot(0:T_last-1,emp_dev_nogrant(1,1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,emp_dev_nogrant(2,1:T_last),'--','LineWidth',LW)
%hold on
%plot(0:T_last-1,emp_dev_nogrant(4,1:T_last),':','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Quartile 1','Quartile 2','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'emp_dev_nogrant'),FormatFig); end

% capital, baseline grant
figure
plot(0:T_last-1,K_dev_grant_baseline(1,1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,K_dev_grant_baseline(2,1:T_last),'--','LineWidth',LW)
%hold on
%plot(0:T_last-1,K_dev_grant_baseline(4,1:T_last),':','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Quartile 1','Quartile 2','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'K_dev_grant_baseline'),FormatFig); end

% capital, no grant
figure
plot(0:T_last-1,K_dev_nogrant(1,1:T_last),'-','LineWidth',LW)
hold on
plot(0:T_last-1,K_dev_nogrant(2,1:T_last),'--','LineWidth',LW)
%hold on
%plot(0:T_last-1,K_dev_nogrant(4,1:T_last),':','LineWidth',LW)
yline(0,'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% change','FontSize',FS)
legend('Quartile 1','Quartile 2','FontSize',FS,'location','best')
set(gca,'FontSize',FS)
if do_save==1; print(fullfile(SaveDir,'K_dev_nogrant'),FormatFig); end

%% Compute fraction of saved firms and zombie firms
clear; close all
disp('Compute fraction of saved firms and zombie firms')
% baseline grant
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'pol_tran','distrib_tran','path','par')
else
    error('File grant_baseline.mat does not exist!')
end

[out_zombie_grant_baseline] = fun_zombie(pol_tran,distrib_tran,path,par);
clearvars -except out_zombie_grant_baseline

% targeted grant (only impacted firms receive grant)
if isfile(fullfile('mat','grant_targslim.mat'))
    load(fullfile('mat','grant_targslim.mat'),'pol_tran','distrib_tran','path','par')
else
    error('File grant_targslim.mat does not exist!')
end
[out_zombie_grant_targslim] = fun_zombie(pol_tran,distrib_tran,path,par);
clearvars -except out_zombie_grant_baseline  out_zombie_grant_targslim

% large targeted grant (only impacted firms receive grant)
if isfile(fullfile('mat','grant_targslim_large.mat'))
    load(fullfile('mat','grant_targslim_large.mat'),'pol_tran','distrib_tran','path','par')
else
    error('File grant_targslim_large.mat does not exist!')
end
[out_zombie_grant_targslim_large] = fun_zombie(pol_tran,distrib_tran,path,par);
clearvars -except out_zombie_grant_baseline  out_zombie_grant_targslim out_zombie_grant_targslim_large

% Display zombie results on screen
disp('--------------------------')
disp('Zombie firms analysis: ')
disp('Impact period only')
disp('----Baseline grant')
fprintf('Fraction of saved firms : %f \n', ...
    (out_zombie_grant_baseline.mass_zombie + out_zombie_grant_baseline.mass_no_zombie)/(out_zombie_grant_baseline.mass_nogrant_exit))
fprintf('Fraction of zombie firms : %f \n', ...
    out_zombie_grant_baseline.mass_zombie/(out_zombie_grant_baseline.mass_zombie+out_zombie_grant_baseline.mass_no_zombie))
disp('----Targeted grant')
fprintf('Fraction of saved firms : %f \n', ...
    (out_zombie_grant_targslim.mass_zombie + out_zombie_grant_targslim.mass_no_zombie)/(out_zombie_grant_targslim.mass_nogrant_exit))
fprintf('Fraction of zombie firms : %f \n', ...
    out_zombie_grant_targslim.mass_zombie/(out_zombie_grant_targslim.mass_zombie+out_zombie_grant_targslim.mass_no_zombie))
disp('----Large Targeted grant')
fprintf('Fraction of saved firms : %f \n', ...
    (out_zombie_grant_targslim_large.mass_zombie + out_zombie_grant_targslim_large.mass_no_zombie)/(out_zombie_grant_targslim_large.mass_nogrant_exit))
fprintf('Fraction of zombie firms : %f \n', ...
    out_zombie_grant_targslim_large.mass_zombie/(out_zombie_grant_targslim_large.mass_zombie+out_zombie_grant_targslim_large.mass_no_zombie))
disp('--------------------------')

% Write to txt file
FID = fopen(fullfile('tables','zombie.txt'),'w');
fprintf(FID,'--------------------------\n');
fprintf(FID,'Zombie firms analysis:  \n');
fprintf(FID,'Impact period only \n');
fprintf(FID,'Baseline grant \n');
fprintf(FID,'Fraction of saved firms: %f \n', ...
    (out_zombie_grant_baseline.mass_zombie + out_zombie_grant_baseline.mass_no_zombie)/(out_zombie_grant_baseline.mass_nogrant_exit));
fprintf(FID,'Fraction of zombie firms: %f \n', ...
    out_zombie_grant_baseline.mass_zombie/(out_zombie_grant_baseline.mass_zombie+out_zombie_grant_baseline.mass_no_zombie));
fprintf(FID,'Targeted grant \n');
fprintf(FID,'Fraction of saved firms: %f \n', ...
    (out_zombie_grant_targslim.mass_zombie + out_zombie_grant_targslim.mass_no_zombie)/(out_zombie_grant_targslim.mass_nogrant_exit));
fprintf(FID,'Fraction of zombie firms: %f \n', ...
    out_zombie_grant_targslim.mass_zombie/(out_zombie_grant_targslim.mass_zombie+out_zombie_grant_targslim.mass_no_zombie));
fprintf(FID,'Large Targeted grant \n');
fprintf(FID,'Fraction of saved firms: %f \n', ...
    (out_zombie_grant_targslim_large.mass_zombie + out_zombie_grant_targslim_large.mass_no_zombie)/(out_zombie_grant_targslim_large.mass_nogrant_exit));
fprintf(FID,'Fraction of zombie firms: %f \n', ...
    out_zombie_grant_targslim_large.mass_zombie/(out_zombie_grant_targslim_large.mass_zombie+out_zombie_grant_targslim_large.mass_no_zombie));
fprintf(FID,'--------------------------');
fclose(FID);

%fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
% Write to LATEX file
FID = fopen(fullfile('tables','zombie.tex'),'w');
fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
%fprintf(FID,'Zombie firms analysis:  \n');
%fprintf(FID,'Impact period only \n');
fprintf(FID,'Baseline grant & \\\\ \n');
fprintf(FID,'Fraction of saved firms & %8.4f \\\\ \n', ...
    (out_zombie_grant_baseline.mass_zombie + out_zombie_grant_baseline.mass_no_zombie)/(out_zombie_grant_baseline.mass_nogrant_exit));
fprintf(FID,'Fraction of zombie firms & %8.4f \\\\ \n', ...
    out_zombie_grant_baseline.mass_zombie/(out_zombie_grant_baseline.mass_zombie+out_zombie_grant_baseline.mass_no_zombie));
fprintf(FID, '\\hline \n');
fprintf(FID,'Targeted grant & \\\\ \n');
fprintf(FID,'Fraction of saved firms & %8.4f \\\\  \n', ...
    (out_zombie_grant_targslim.mass_zombie + out_zombie_grant_targslim.mass_no_zombie)/(out_zombie_grant_targslim.mass_nogrant_exit));
fprintf(FID,'Fraction of zombie firms & %8.4f \\\\ \n', ...
    out_zombie_grant_targslim.mass_zombie/(out_zombie_grant_targslim.mass_zombie+out_zombie_grant_targslim.mass_no_zombie));
fprintf(FID,'Targeted grant (large) & \\\\ \n');
fprintf(FID,'Fraction of saved firms & %8.4f \\\\  \n', ...
    (out_zombie_grant_targslim_large.mass_zombie + out_zombie_grant_targslim_large.mass_no_zombie)/(out_zombie_grant_targslim_large.mass_nogrant_exit));
fprintf(FID,'Fraction of zombie firms & %8.4f \\\\ \n', ...
    out_zombie_grant_targslim_large.mass_zombie/(out_zombie_grant_targslim_large.mass_zombie+out_zombie_grant_targslim_large.mass_no_zombie));
fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);

%% Four types of capital adjustment
disp('Four types of capital adjustment')
clear; close all
% Load ss results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'agg','par')
else
    error('File ss.mat does not exist!')
end

capadj_ss = agg.capadj;
% Flip the sign to negative for downward adjustment and capital sold by
% exiters.
capadj_ss([2, 4]) = - capadj_ss([2, 4]);
% Capital adjustment rate
capadj_rate_ss = capadj_ss/(agg.K_small);

% Write to LATEX file
FID = fopen(fullfile('tables','capadj_rate_ss.tex'),'w');
fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
%fprintf(FID,'Zombie firms analysis:  \n');
%fprintf(FID,'Impact period only \n');
fprintf(FID,'Up. adj. by active firms & %8.4f \\\\ \n',capadj_rate_ss(1));
fprintf(FID,'Down. adj. by active firms & %8.4f \\\\ \n',capadj_rate_ss(2));
fprintf(FID,'Capital bought by entrants & %8.4f \\\\ \n',capadj_rate_ss(3));
fprintf(FID,'Capital sold by exiters & %8.4f \\\\ \n',capadj_rate_ss(4));
fprintf(FID,'Overall small-firm investment rate & %8.4f \\\\ \n',...
    sum(capadj_rate_ss));
fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);


% Load transition results
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'agg_tran','irf')
else
    error('File grant_baseline.mat does not exist!')
end

capadj_baseline = agg_tran.capadj;
% Flip the sign to negative for downward adjustment and capital sold by
% exiters.
capadj_baseline(:,[2,4]) = -capadj_baseline(:,[2,4]);
% Continuing firms' investment decisions are only realized in the
% subsequent period.
capadj_baseline(2:par.T+1,1:2) =  capadj_baseline(1:par.T,1:2);
capadj_baseline(1,1:2) = capadj_ss(1:2)';
capadj_rate_baseline = capadj_baseline/agg.K_small;
%capadj_baseline_diff = 100*(capadj_baseline - capadj_ss')./agg.K_small;
%capadj_share_baseline = capadj_baseline./abs(capadj_baseline(:,1)-capadj_baseline(:,2)+capadj_baseline(:,3)-capadj_baseline(:,4));

if isfile(fullfile('mat','nogrant.mat'))
    load(fullfile('mat','nogrant.mat'),'agg_tran')
else
    error('File nogrant.mat does not exist!')
end
capadj_nogrant = agg_tran.capadj;
% Flip the sign to negative for downward adjustment and capital sold by
% exiters.
capadj_nogrant(:,[2,4]) = -capadj_nogrant(:,[2,4]);
% Continuing firms' investment decisions are only realized in the
% subsequent period.
capadj_nogrant(2:par.T+1,1:2) =  capadj_nogrant(1:par.T,1:2);
capadj_nogrant(1,1:2) = capadj_ss(1:2);
capadj_rate_nogrant = capadj_nogrant/agg.K_small;
%capadj_nogrant_diff = 100*(capadj_nogrant - capadj_ss')./capadj_ss';

if isfile(fullfile('mat','grant_targslim.mat'))
    load(fullfile('mat','grant_targslim.mat'),'agg_tran')
else
    error('File grant_targslim.mat does not exist!')
end
capadj_targslim = agg_tran.capadj;
% Flip the sign to negative for downward adjustment and capital sold by
% exiters.
capadj_targslim(:,[2,4]) = -capadj_targslim(:,[2,4]);
% Continuing firms' investment decisions are only realized in the
% subsequent period.
capadj_targslim(2:par.T+1,1:2) =  capadj_targslim(1:par.T,1:2);
capadj_targslim(1,1:2) = capadj_ss(1:2);
capadj_rate_targslim = capadj_targslim/agg.K_small;
%capadj_targslim_diff = 100*(capadj_targslim - capadj_ss')./capadj_ss';

if isfile(fullfile('mat','grant_targslim_large.mat'))
    load(fullfile('mat','grant_targslim_large.mat'),'agg_tran')
else
    error('File grant_targslim_large.mat does not exist!')
end
capadj_targslim_large = agg_tran.capadj;
% Flip the sign to negative for downward adjustment and capital sold by
% exiters.
capadj_targslim_large(:,[2,4]) = -capadj_targslim_large(:,[2,4]);
% Continuing firms' investment decisions are only realized in the
% subsequent period.
capadj_targslim_large(2:par.T+1,1:2) =  capadj_targslim_large(1:par.T,1:2);
capadj_targslim_large(1,1:2) = capadj_ss(1:2);
capadj_rate_targslim_large = capadj_targslim_large/agg.K_small;
%capadj_targslim_diff = 100*(capadj_targslim - capadj_ss')./capadj_ss';


% Number of quarters in the short-run, mid-run, and long-run
T_sr = 2; %
T_mr = 8; %
T_lr = 40; % 10 years
% Create a stacked bar chart using the bar function
sr_cum_capadj_nogrant         = zeros(4,1);
sr_cum_capadj_grant_baseline  = zeros(4,1);
sr_cum_capadj_grant_targslim  = zeros(4,1);
sr_cum_capadj_grant_targslim_large  = zeros(4,1);
mr_cum_capadj_nogrant         = zeros(4,1);
mr_cum_capadj_grant_baseline  = zeros(4,1);
mr_cum_capadj_grant_targslim  = zeros(4,1);
mr_cum_capadj_grant_targslim_large  = zeros(4,1);
lr_cum_capadj_nogrant         = zeros(4,1);
lr_cum_capadj_grant_baseline  = zeros(4,1);
lr_cum_capadj_grant_targslim  = zeros(4,1);
lr_cum_capadj_grant_targslim_large  = zeros(4,1);
ave_cum_capadj_nogrant        = zeros(4,1);
ave_cum_capadj_grant_baseline = zeros(4,1);
ave_cum_capadj_grant_targslim = zeros(4,1);
ave_cum_capadj_grant_targslim_large = zeros(4,1);

for ii = 1:4

    sr_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(1:T_sr,ii)-capadj_rate_ss(ii));
    sr_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(1:T_sr,ii)-capadj_rate_ss(ii));
    sr_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(1:T_sr,ii)-capadj_rate_ss(ii));
    sr_cum_capadj_grant_targslim_large(ii) = 100*sum(capadj_rate_targslim_large(1:T_sr,ii)-capadj_rate_ss(ii));

    mr_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(T_sr+1:T_mr,ii)-capadj_rate_ss(ii));
    mr_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(T_sr+1:T_mr,ii)-capadj_rate_ss(ii));
    mr_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(T_sr+1:T_mr,ii)-capadj_rate_ss(ii));
    mr_cum_capadj_grant_targslim_large(ii) = 100*sum(capadj_rate_targslim_large(T_sr+1:T_mr,ii)-capadj_rate_ss(ii));

    lr_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(T_mr+1:T_lr,ii)-capadj_rate_ss(ii));
    lr_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(T_mr+1:T_lr,ii)-capadj_rate_ss(ii));
    lr_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(T_mr+1:T_lr,ii)-capadj_rate_ss(ii));
    lr_cum_capadj_grant_targslim_large(ii) = 100*sum(capadj_rate_targslim_large(T_mr+1:T_lr,ii)-capadj_rate_ss(ii));

    ave_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(1:T_lr,ii)-capadj_rate_ss(ii));
    ave_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(1:T_lr,ii)-capadj_rate_ss(ii));
    ave_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(1:T_lr,ii)-capadj_rate_ss(ii));
    ave_cum_capadj_grant_targslim_large(ii) = 100*sum(capadj_rate_targslim_large(1:T_lr,ii)-capadj_rate_ss(ii));

end
% Get the sum over the 4 types of capital adjustment, for each policy
% environ. and for each period (short run, etc)
sum_sr_cum_capadj_nogrant  = sum(sr_cum_capadj_nogrant);
sum_mr_cum_capadj_nogrant  = sum(mr_cum_capadj_nogrant);
sum_lr_cum_capadj_nogrant  = sum(lr_cum_capadj_nogrant);
sum_ave_cum_capadj_nogrant = sum(ave_cum_capadj_nogrant);

sum_sr_cum_capadj_grant_baseline  = sum(sr_cum_capadj_grant_baseline);
sum_mr_cum_capadj_grant_baseline  = sum(mr_cum_capadj_grant_baseline);
sum_lr_cum_capadj_grant_baseline  = sum(lr_cum_capadj_grant_baseline);
sum_ave_cum_capadj_grant_baseline = sum(ave_cum_capadj_grant_baseline);

sum_sr_cum_capadj_grant_targslim  = sum(sr_cum_capadj_grant_targslim);
sum_mr_cum_capadj_grant_targslim  = sum(mr_cum_capadj_grant_targslim);
sum_lr_cum_capadj_grant_targslim  = sum(lr_cum_capadj_grant_targslim);
sum_ave_cum_capadj_grant_targslim = sum(ave_cum_capadj_grant_targslim);

sum_sr_cum_capadj_grant_targslim_large  = sum(sr_cum_capadj_grant_targslim_large);
sum_mr_cum_capadj_grant_targslim_large  = sum(mr_cum_capadj_grant_targslim_large);
sum_lr_cum_capadj_grant_targslim_large  = sum(lr_cum_capadj_grant_targslim_large);
sum_ave_cum_capadj_grant_targslim_large = sum(ave_cum_capadj_grant_targslim_large);


% make plots
T_last    = 40;
FS        = 18;
LW        = 3;
do_save   = 1;
FormatFig = '-depsc';  % Specify '-dpng' or '-depsc'
SaveDir   = fullfile('figures','grant_vs_nogrant');


fileroot = '';
plot_capadj(capadj_rate_ss,capadj_rate_baseline,capadj_rate_nogrant,capadj_rate_targslim,...
    sr_cum_capadj_nogrant,sr_cum_capadj_grant_baseline,sr_cum_capadj_grant_targslim,...
    sum_sr_cum_capadj_nogrant,sum_sr_cum_capadj_grant_baseline,sum_sr_cum_capadj_grant_targslim,...
    mr_cum_capadj_nogrant,mr_cum_capadj_grant_baseline,mr_cum_capadj_grant_targslim,...
    sum_mr_cum_capadj_nogrant,sum_mr_cum_capadj_grant_baseline,sum_mr_cum_capadj_grant_targslim,...
    lr_cum_capadj_nogrant,lr_cum_capadj_grant_baseline,lr_cum_capadj_grant_targslim,...
    sum_lr_cum_capadj_nogrant,sum_lr_cum_capadj_grant_baseline,sum_lr_cum_capadj_grant_targslim,...
    ave_cum_capadj_nogrant,ave_cum_capadj_grant_baseline,ave_cum_capadj_grant_targslim,...
    sum_ave_cum_capadj_nogrant,sum_ave_cum_capadj_grant_baseline,sum_ave_cum_capadj_grant_targslim,...
    T_last,LW,FS,SaveDir,fileroot,FormatFig,do_save)
fileroot = '_large';
plot_capadj(capadj_rate_ss,capadj_rate_baseline,capadj_rate_nogrant,capadj_rate_targslim_large,...
    sr_cum_capadj_nogrant,sr_cum_capadj_grant_baseline,sr_cum_capadj_grant_targslim_large,...
    sum_sr_cum_capadj_nogrant,sum_sr_cum_capadj_grant_baseline,sum_sr_cum_capadj_grant_targslim_large,...
    mr_cum_capadj_nogrant,mr_cum_capadj_grant_baseline,mr_cum_capadj_grant_targslim_large,...
    sum_mr_cum_capadj_nogrant,sum_mr_cum_capadj_grant_baseline,sum_mr_cum_capadj_grant_targslim_large,...
    lr_cum_capadj_nogrant,lr_cum_capadj_grant_baseline,lr_cum_capadj_grant_targslim_large,...
    sum_lr_cum_capadj_nogrant,sum_lr_cum_capadj_grant_baseline,sum_lr_cum_capadj_grant_targslim_large,...
    ave_cum_capadj_nogrant,ave_cum_capadj_grant_baseline,ave_cum_capadj_grant_targslim_large,...
    sum_ave_cum_capadj_nogrant,sum_ave_cum_capadj_grant_baseline,sum_ave_cum_capadj_grant_targslim_large,...
    T_last,LW,FS,SaveDir,fileroot,FormatFig,do_save)

%% Welfare analysis
clear
% Load ss results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'agg')
else
    error('File ss.mat does not exist!')
end

C_ss = agg.C_agg;
L_ss = agg.L_agg;

% Load transition results
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'agg_tran','path','par')
else
    error('File grant_baseline.mat does not exist!')
end
C_baseline = path.C;
L_baseline = agg_tran.L_agg;
margutil   = par.margutil;
lsupply    = par.lsupply;

if isfile(fullfile('mat','nogrant.mat'))
    load(fullfile('mat','nogrant.mat'),'agg_tran','path')
else
    error('File nogrant.mat does not exist!')
end
C_nogrant = path.C;
L_nogrant = agg_tran.L_agg;

if isfile(fullfile('mat','grant_targslim.mat'))
    load(fullfile('mat','grant_targslim.mat'),'agg_tran','path')
else
    error('File grant_targslim.mat does not exist!')
end
C_targslim = path.C;
L_targslim = agg_tran.L_agg;

if isfile(fullfile('mat','grant_targslim_large.mat'))
    load(fullfile('mat','grant_targslim_large.mat'),'agg_tran','path')
else
    error('File grant_targslim_large.mat does not exist!')
end
C_targslim_large = path.C;
L_targslim_large = agg_tran.L_agg;

[CEV_baseline,CEV_targslim] = fun_welfare(par,C_ss,L_ss,C_baseline,L_baseline,...
    C_nogrant,L_nogrant,C_targslim,L_targslim,margutil,lsupply);
[~,CEV_targslim_large] = fun_welfare(par,C_ss,L_ss,C_baseline,L_baseline,...
    C_nogrant,L_nogrant,C_targslim_large,L_targslim_large,margutil,lsupply);

% CEV in terms of steady state GDP(?)
%CEV_baseline*C_baseline(1)/(1-par.beta)/agg.Y_agg


% Write to txt file
FID = fopen(fullfile('tables','welfare_analysis.txt'),'w');
fprintf(FID,'--------------------------\n');
fprintf(FID,'CEV: grant vs. no grant \n');
fprintf(FID,'Baseline grant: %f \n',CEV_baseline);
fprintf(FID,'Targeted grant: %f \n', CEV_targslim);
fprintf(FID,'Large Targeted grant: %f \n', CEV_targslim_large);
fprintf(FID,'--------------------------\n');
fclose(FID);

% Write to Latex file
FID = fopen(fullfile('tables','welfare_analysis.tex'),'w');
fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
%fprintf(FID,'CEV: grant vs. no grant &       \\\\ \n');
fprintf(FID,'Baseline grant         & %8.6f \\\\ \n', CEV_baseline);
fprintf(FID,'Targeted grant         & %8.6f \\\\ \n', CEV_targslim);
fprintf(FID,'Targeted grant (large) & %8.6f \\\\ \n', CEV_targslim_large);
fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);

