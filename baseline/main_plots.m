%% This script generates figures (and some tables/other results) for the paper
% IMPORTANT: Before running this file you must have saved results of the
% model (run main.m before you run this file).
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

% Set up some useful paths
FormatFig = '-dpng';  % Specify '-dpng' or '-depsc'

matNames = {'ss','nogrant','grant_baseline','grant_targslim'};
for ii = 1:numel(matNames)
    if ~isfile(fullfile('mat',[matNames{ii},'.mat']))
        warning('MAT file ''%s'' is missing \n',matNames{ii})
    end
end

%% Plot Steady-state distributions, policy functions (exit, entry, investment)
disp('Plot Steady-state distributions, policy functions')
% Load results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'b_grid','sol','par','prices','distribS','model_mom','data_mom') %#ok<NASGU> 
else
    error('File ss.mat does not exist!')
end

% Specify here where you want to save the figures
SaveDir = fullfile('figures','ss');

FS  = 16; % font size for plots

% Call function to make plots and to save them
plot_ss(b_grid,sol,par,distribS,model_mom,data_mom,SaveDir,1,FormatFig,FS);

% We plot policy function of investment
plot_ss_policy(b_grid,sol,distribS,par,SaveDir,FormatFig,1,FS)

%% Plot calibration of transition, and plot the exit thresholds in t=1 for the baseline grant
clear; close all
disp('Plot calibration of transition')
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'model_mom_trans','par');%,'b_grid','sol','pol_tran');
else
    error('File grant_baseline.mat does not exist')
end

FormatFig = '-dpng';  % Specify '-dpng' or '-depsc'
FS        = 16; % font size for plots
do_save   = 1;
SaveDir   = fullfile('figures','grant_vs_nogrant');
% Calibration of transition: moments and shocks
plot_calib_tran(model_mom_trans,par,FormatFig,FS,do_save,SaveDir);

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
FormatFig = '-dpng'; % Specify '-dpng' or '-depsc'
FS = 16; % font size for plots
LW = 3; % line width
T_last = min(70,par.T);
do_save = 1;
SaveDir = fullfile('figures','grant_vs_nogrant');

% Do plots
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

nbins = 4;
[exitSS,exit_nogrant,exit_grant_baseline,binlids] = plot_micro_l(sol,distribS,...
    prices,pol_tran_nogrant,weights_nogrant,pol_tran_grant_baseline,...
    weights_grant_baseline,par);

% Change in exit rate from no grant economy to grant economy
change_exit_rate  = exit_grant_baseline.exit_rate - exit_nogrant.exit_rate;
change_exit       = exit_grant_baseline.exits - exit_nogrant.exits;
share_change_exit = change_exit / sum(change_exit);

% Options for plots
FS        = 13; % font size
FORMAT    = '-dpng'; % Specify '-dpng' or '-depsc'
SaveDir   = fullfile('figures','grant_vs_nogrant');
do_save   = 1;
bin_names = {'0-9';'10-19';'20-99';'100+'};

% Employment share in ss vs share of change in exit due to grant, by decile of x
figure
yyaxis left
w1 = 0.8;
bar(1:nbins,exitSS.emp_share,w1)
hold on
yyaxis right
w2 = 0.6;
bar(1:nbins,change_exit_rate,w2)
set(gca,'FontSize',12);
set(gca,'xticklabel',bin_names)
xlabel('Bins by firm size (employment)','FontSize',FS)
%title('Effect of grant on exit rate','FontSize',14)
legend('Employment share in steady state (left axis)','Change in exit rate (right axis)','FontSize',FS,'location','southoutside')
hold off
if do_save==1; print(fullfile(SaveDir,'emp_exit_rate_lbins'),FORMAT); end

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

% Write results to LATEX file
FID = fopen(fullfile('tables','zombie.tex'),'w');
fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
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
fprintf(FID,'Up. adj. by active firms & %8.3f \\\\ \n',capadj_rate_ss(1));
fprintf(FID,'Down. adj. by active firms & %8.3f \\\\ \n',capadj_rate_ss(2));
fprintf(FID,'Capital bought by entrants & %8.3f \\\\ \n',capadj_rate_ss(3));
fprintf(FID,'Capital sold by exiters & %8.3f \\\\ \n',capadj_rate_ss(4));
fprintf(FID,'Overall small-firm investment rate & %8.3f \\\\ \n',...
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

% Make bar plots
T_last    = 40;
FS        = 20;
LW        = 3;
do_save   = 1;
FormatFig = '-dpng';  % Specify '-dpng' or '-depsc'
SaveDir   = fullfile('figures','grant_vs_nogrant');

% Number of quarters in the short-run, mid-run, and long-run
T_sr = 2; %
T_mr = 8; %
T_lr = 40; % 10 years
% Create a stacked bar chart using the bar function
sr_cum_capadj_nogrant         = zeros(4,1);
sr_cum_capadj_grant_baseline  = zeros(4,1);
sr_cum_capadj_grant_targslim  = zeros(4,1);
mr_cum_capadj_nogrant         = zeros(4,1);
mr_cum_capadj_grant_baseline  = zeros(4,1);
mr_cum_capadj_grant_targslim  = zeros(4,1);
lr_cum_capadj_nogrant         = zeros(4,1);
lr_cum_capadj_grant_baseline  = zeros(4,1);
lr_cum_capadj_grant_targslim  = zeros(4,1);
ave_cum_capadj_nogrant        = zeros(4,1);
ave_cum_capadj_grant_baseline = zeros(4,1);
ave_cum_capadj_grant_targslim = zeros(4,1);

for ii = 1:4

    sr_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(1:T_sr,ii)-capadj_rate_ss(ii));
    sr_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(1:T_sr,ii)-capadj_rate_ss(ii));
    sr_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(1:T_sr,ii)-capadj_rate_ss(ii));

    mr_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(T_sr+1:T_mr,ii)-capadj_rate_ss(ii));
    mr_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(T_sr+1:T_mr,ii)-capadj_rate_ss(ii));
    mr_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(T_sr+1:T_mr,ii)-capadj_rate_ss(ii));

    lr_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(T_mr+1:T_lr,ii)-capadj_rate_ss(ii));
    lr_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(T_mr+1:T_lr,ii)-capadj_rate_ss(ii));
    lr_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(T_mr+1:T_lr,ii)-capadj_rate_ss(ii));

    ave_cum_capadj_nogrant(ii)        = 100*sum(capadj_rate_nogrant(1:T_lr,ii)-capadj_rate_ss(ii));
    ave_cum_capadj_grant_baseline(ii) = 100*sum(capadj_rate_baseline(1:T_lr,ii)-capadj_rate_ss(ii));
    ave_cum_capadj_grant_targslim(ii) = 100*sum(capadj_rate_targslim(1:T_lr,ii)-capadj_rate_ss(ii));

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

% Cumulative capital adjstment in the short run: t = 1 to T_sr
figure('Position',[1000 918 560 420])
b=bar(1:3, [sr_cum_capadj_nogrant, sr_cum_capadj_grant_baseline, sr_cum_capadj_grant_targslim]', 0.5, 'stack'); %#ok<NASGU> 
hold on
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Baseline grant' ...
    'Targeted grant' },'fontsize',FS)
xtickangle(30)
hL=plot(1:3,[sum_sr_cum_capadj_nogrant sum_sr_cum_capadj_grant_baseline ...
    sum_sr_cum_capadj_grant_targslim],'sg','MarkerSize',15,'MarkerFaceColor','g'); %#ok<NASGU> 
hold off
%title(varlabel,'fontsize',myfontsize)
ylabel({'Cumulative excess investment', '(% of steady state capital)'},'fontsize',FS)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = FS;
legend('Invest', 'Disinvest','Entry', 'Exit','All adj.','location','northoutside',...
    'fontsize',FS-1,'NumColumns',5)
print(fullfile(SaveDir,'cum_capadj_decomp_sr'),FormatFig);

% Cumulative capital adjstment in the medium run: t = T_sr+1 to T_mr
figure('Position',[1000 918 560 420])
b=bar(1:3, [mr_cum_capadj_nogrant, mr_cum_capadj_grant_baseline, mr_cum_capadj_grant_targslim]', 0.5, 'stack'); %#ok<NASGU> 
hold on
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Baseline grant' ...
    'Targeted grant' },'fontsize',FS)
xtickangle(30)
hL=plot(1:3,[sum_mr_cum_capadj_nogrant sum_mr_cum_capadj_grant_baseline ...
    sum_mr_cum_capadj_grant_targslim],'sg','MarkerSize',15,'MarkerFaceColor','g'); %#ok<NASGU> 
hold off
%title(varlabel,'fontsize',myfontsize)
ylabel({'Cumulative excess investment', '(% of steady state capital)'},'fontsize',FS)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = FS;
legend('Invest', 'Disinvest','Entry', 'Exit','All adj.','location','northoutside',...
    'fontsize',FS-1,'NumColumns',5)
print(fullfile(SaveDir,'cum_capadj_decomp_mr'),FormatFig);

% Cumulative capital adjstment in the long run: t = T_mr+1 to T_lr
figure('Position',[1000 918 560 420])
b=bar(1:3, [lr_cum_capadj_nogrant, lr_cum_capadj_grant_baseline, lr_cum_capadj_grant_targslim]', 0.5, 'stack'); %#ok<NASGU> 
hold on
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Baseline grant' ...
    'Targeted grant' },'fontsize',FS)
xtickangle(30)
hL=plot(1:3,[sum_lr_cum_capadj_nogrant sum_lr_cum_capadj_grant_baseline ...
    sum_lr_cum_capadj_grant_targslim],'sg','MarkerSize',15,'MarkerFaceColor','g'); %#ok<NASGU> 
hold off
%title(varlabel,'fontsize',myfontsize)
ylabel({'Cumulative excess investment', '(% of steady state capital in small firms)'},'fontsize',FS)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = FS;
legend('Invest', 'Disinvest','Entry', 'Exit','All adj.','location','northoutside',...
    'fontsize',FS-1,'NumColumns',5)
print(fullfile(SaveDir,'cum_capadj_decomp_lr'),FormatFig);

% Cumulative capital adjstment: t = 1 to T_lr
figure('Position',[1000 918 560 420])
b=bar(1:3, [ave_cum_capadj_nogrant, ave_cum_capadj_grant_baseline, ave_cum_capadj_grant_targslim]', 0.5, 'stack');
hold on
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Baseline grant' ...
    'Targeted grant' },'fontsize',FS)
xtickangle(30)
hL=plot(1:3,[sum_ave_cum_capadj_nogrant sum_ave_cum_capadj_grant_baseline ...
    sum_ave_cum_capadj_grant_targslim],'sg','MarkerSize',15,'MarkerFaceColor','g');
hold off
%title(varlabel,'fontsize',myfontsize)
ylabel({'Cumulative excess investment', '(% of steady state capital)'},'fontsize',FS)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = FS;
legend('Invest', 'Disinvest','Entry', 'Exit','All adj.','location','northoutside',...
    'fontsize',FS-1,'NumColumns',5)
print(fullfile(SaveDir,'cum_capadj_decomp_all'),FormatFig);

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

[CEV_baseline,CEV_targslim] = fun_welfare(par,C_ss,L_ss,C_baseline,L_baseline,...
    C_nogrant,L_nogrant,C_targslim,L_targslim,margutil,lsupply);

% Write to txt file
FID = fopen(fullfile('tables','welfare_analysis.txt'),'w');
fprintf(FID,'--------------------------\n');
fprintf(FID,'CEV: grant vs. no grant \n');
fprintf(FID,'Baseline grant: %f \n',CEV_baseline);
fprintf(FID,'Targeted grant: %f \n', CEV_targslim);
fprintf(FID,'--------------------------\n');
fclose(FID);

% Write to Latex file
FID = fopen(fullfile('tables','welfare_analysis.tex'),'w');
fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
%fprintf(FID,'CEV: grant vs. no grant &       \\\\ \n');
fprintf(FID,'Baseline grant        & %8.6f \\\\ \n', CEV_baseline);
fprintf(FID,'Targeted grant         & %8.6f \\\\ \n', CEV_targslim);
fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);
