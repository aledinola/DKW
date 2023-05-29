%% Make plots to compare large vs small grants. These are uniform grants as in the baseline grant case.
% loading saved results in mat files from subfolder 'mat'

clear;clc;close all

disp('Make plots to compare large vs small grant economy')

%% Set up some useful paths
ResultsDir = fullfile('mat'); %folder where .mat files are stored
FormatFig = '-dpng';  % Specify '-dpng' or '-depsc'


%% Make IRFs for no grant economy
% Load results for no grant economy
load(fullfile(ResultsDir,'nogrant.mat'))

disp('NO GRANT ECONOMY:')

% Compute IRFs in percentage deviation from steady-state
[~,irf_nogrant] = fun_targets_tran(data_mom_trans,agg_tran,path,agg,prices,calibWeightsTran);

clearvars -except ResultsDir FormatFig irf_nogrant

%% Make IRFs for the large grant economy
load(fullfile(ResultsDir,'grant_large.mat'))

disp('LARGE GRANT ECONOMY:')

% Compute IRFs in percentage deviation from steady-state
[~,irf_grant_large] = fun_targets_tran(data_mom_trans,agg_tran,path,agg,prices,calibWeightsTran);

%% Make IRFs for small grant economy
load(fullfile(ResultsDir,'grant_small.mat'))

disp('SMALL GRANT ECONOMY:')

% Compute IRFs in percentage deviation from steady-state
[~,irf_grant_small] = fun_targets_tran(data_mom_trans,agg_tran,path,agg,prices,calibWeightsTran);


%% Plot IRFs of economy comparing no grant, baseline grant and targeted grant
% Specify here where you want to save the figures

% Call function to make plots and to save them
% plot_irf_compare(irf1,irf2,irf_nogrant,name1,name2,data_mom_trans,SaveDir,do_save,FormatFig)
name1 = 'Large grant';
name2 = 'Small grant';
% Specify here where you want to save the figures
FigDir = fullfile('figures','large_vs_small');
if ~(isfolder(FigDir))
    disp("Folder does not exist, creating it now..")
    mkdir(FigDir)
end
plot_irf_compare(irf_grant_large,irf_grant_small,irf_nogrant,name1,name2,data_mom_trans,FigDir,1,FormatFig);


