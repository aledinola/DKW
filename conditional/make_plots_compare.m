%% Make plots to compare no grant vs grant economy
% We need the following mat files, stored in folder mat:
% - ss.mat
% - nogrant.mat
% - grant_baseline.mat
% - grant_targslim.mat
% This script calls the following functions:
% - plot_irf_compare.m

clear;clc;close all

matNames = {'ss','nogrant','grant_baseline'};
for ii = 1:numel(matNames)
    if ~isfile(fullfile('mat',[matNames{ii},'.mat']))
        warning('MAT file ''%s'' is missing \n',matNames{ii})
    end
end


disp('Make plots to compare no grant vs grant economy')

%% Set up some useful paths
ResultsDir = fullfile('mat'); %folder where .mat files are stored
% Specify here where you want to save the figures

%% Load ss results
% Load results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'agg')
else
    error('File ss.mat does not exist!')
end
agg_ss = agg;
clear agg

%% Load IRFs 
% Load results for no grant economy
load(fullfile(ResultsDir,'nogrant.mat'),'irf','agg_tran')

disp('NO GRANT ECONOMY:')
irf_nogrant = irf;
ppchange_nogrant.exit_rate = agg_tran.exit_rate_vec - agg_ss.exit_rate;
ppchange_nogrant.entry_rate = agg_tran.entry_rate_vec - agg_ss.entry_rate;

% Compute IRFs in percentage deviation from steady-state
%[~,irf_nogrant] = fun_targets_tran(data_mom_trans,agg_tran,path,agg,prices,calibWeightsTran);
% Compute IRFs by impact status
%[irf_impact_nogrant] = fun_aggregates_tran_impact(par,pol_tran,distrib_tran,path,agg);

clear irf agg_tran

%% Make IRFs for baseline grant economy
load(fullfile(ResultsDir,'grant_baseline.mat'),'irf','agg_tran')

disp('BASELINE GRANT ECONOMY:')
irf_grant_baseline = irf;
ppchange_grant_baseline.exit_rate = agg_tran.exit_rate_vec - agg_ss.exit_rate;
ppchange_grant_baseline.entry_rate = agg_tran.entry_rate_vec - agg_ss.entry_rate;

clear irf agg_tran

%% Plot IRFs 

% Options for plots
FormatFig = '-dpng'; % Specify '-dpng' or '-depsc'
do_save   = 1;        % Flag 0-1
lastp     = 40;       % number of periods plotted
% - irf: Baseline grant vs laissez-faire
name1 = 'Conditional grant';%'Baseline grant';
name2 = 'Laissez-faire';
plot_irf_compare(irf_grant_baseline,irf_nogrant,...
    name1,name2,lastp,fullfile('figures','grant_vs_nogrant'),do_save,FormatFig);

% percentage point change of entry and exit rates, baseline grant vs laissez-faire

plot_ppchange_compare(ppchange_grant_baseline,ppchange_nogrant,...
    name1,name2,lastp,fullfile('figures','grant_vs_nogrant'),do_save,FormatFig);

