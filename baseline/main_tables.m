%% This script compiles results (tables) for the paper
% This script requires the following mat files:
%   - ss.mat, nogrant.mat, grant_baseline.mat, grant_targslim.mat 
% This script creates tables used in the writeup.
% loading saved results in mat files
% Tables created by this file:
%  STEADY-STATE
%   - comp_parameters.tex: Table with Computational parameters
%   - exo_parameters.tex: Table with exogenous parameters
%   - steady_state.tex: Table with some steady-state values
%   - moments.tex: Table with model fit
%   - parameters.tex: Table with estimated parameters
%  TRANSITION
%   - tran_shocks.tex: Table with shocks
%   - transition_moments.tex: table with transition moments
%   - cost_grants.tex: table with cost per job saved
%   - cost_grants_robust.tex

clear;clc;close all

matNames = {'ss','nogrant','grant_baseline','grant_targslim'};
for ii = 1:numel(matNames)
    if ~isfile(fullfile('mat',[matNames{ii},'.mat']))
        warning('MAT file ''%s'' is missing \n',matNames{ii})
    end
end


%% Make tables for steady-state

% Load steady-state results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'))
else
    error('File ss.mat does not exist!')
end

% Table with Computational parameters
FID = fopen(fullfile('tables','comp_parameters.tex'),'w');
fprintf(FID,' \\begin{tabular}{llc} \\hline \n');
fprintf(FID,' Parameter & Description & Value \\\\ \n');
fprintf(FID,' \\hline \n');
fprintf(FID,'nx  & Num. of grid points for $x$       & %d \\\\ \n',par.nx);
fprintf(FID,'nb  & Num. of grid points for $b$       & %d \\\\ \n',par.nb);
fprintf(FID,'nk  & Num. of grid points for $\\kappa$  & %d \\\\ \n',par.nk);
fprintf(FID,'T   & Length of transition              & %d \\\\ \n',par.T);
fprintf(FID, '\\hline \n \\end{tabular} \n');
fclose(FID);

% Table with exogenous parameters
FID = fopen(fullfile('tables','exo_parameters.tex'),'w');
fprintf(FID,' \\begin{tabular}{llc} \\hline \n');
fprintf(FID,' Parameter & Description & Value \\\\ \n');
fprintf(FID,' \\hline \n');
for i = 1:numel(ExoNames(:,1))
    fprintf(FID,'%s  & %s & %8.3f \\\\ \n',ExoNames{i,2},ExoNames{i,3},par.(ExoNames{i,1}));
end
fprintf(FID,'%s  & %s & %8.3f \\\\ \n','$bk0(1)$','Initial debt-asset ratio',par.bk0_vec(1));
fprintf(FID,'%s  & %s & %8.3f \\\\ \n','$bk0(2)$','Initial debt-asset ratio',par.bk0_vec(2));
fprintf(FID,'%s  & %s & %8.3f \\\\ \n','$bk0(3)$','Initial debt-asset ratio',par.bk0_vec(3));
fprintf(FID,'%s  & %s & %8.3f \\\\ \n','$Pr(bk0(1))$','Prob. dist.',par.bk0_prob(1));
fprintf(FID,'%s  & %s & %8.3f \\\\ \n','$Pr(bk0(2))$','Prob. dist.',par.bk0_prob(2));
fprintf(FID,'%s  & %s & %8.3f \\\\ \n','$Pr(bk0(3))$','Prob. dist.',par.bk0_prob(3));

fprintf(FID, '\\hline \n \\end{tabular} \n');
fclose(FID);

% Table with estimated parameters (the same for both economies)
% mystruct2table(pstruct,pnames,description,dispNames,header,tex,tabDir,filename)
mystruct2table(par,calibNames,description,dispNames,{'Parameter','Value'},1,fullfile('tables'),'parameters.tex');

% Table with some steady-state values
make_table_ss(prices,agg,b_grid,1,fullfile('tables'),'steady_state.tex');

% Table with model fit
mystruct2table_mom(data_mom,model_mom,targetNames,calibWeights,targetNames_long,1,fullfile('tables'),'moments.tex');

%keyboard

%% Make table for transition
clear

% Load results grant baseline
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'irf','model_mom_trans')
else
    error('File grant_baseline.mat does not exist!')
end
irf_grant_baseline = irf;
model_mom_trans_baseline = model_mom_trans;
clear irf model_mom_trans


% Load results no grant economy
if isfile(fullfile('mat','nogrant.mat'))
    load(fullfile('mat','nogrant.mat'),'par','irf','model_mom_trans','data_mom_trans')
else
    error('File nogrant.mat does not exist!')
end
irf_nogrant = irf;
model_mom_trans_nogrant = model_mom_trans;
clear irf model_mom_trans

% Annual change in employment and output
T_last = 40;%par.T;
tot_Y_agg_grant_baseline = 100*sum(irf_grant_baseline.Y_agg(1:T_last))/T_last;
tot_C_agg_grant_baseline = 100*sum(irf_grant_baseline.C_agg(1:T_last))/T_last;
tot_L_agg_grant_baseline = 100*sum(irf_grant_baseline.L_agg(1:T_last))/T_last;
tot_Y_corp_grant_baseline = 100*sum(irf_grant_baseline.Y_corp(1:T_last))/T_last;
tot_L_corp_grant_baseline = 100*sum(irf_grant_baseline.L_corp(1:T_last))/T_last;

tot_L_small_grant_baseline = 100*sum(irf_grant_baseline.L_small(1:T_last))/T_last;
tot_output_small_grant_baseline = 100*sum(irf_grant_baseline.output_small(1:T_last))/T_last;

tot_Y_agg_nogrant = 100*sum(irf_nogrant.Y_agg(1:T_last))/T_last;
tot_C_agg_nogrant = 100*sum(irf_nogrant.C_agg(1:T_last))/T_last;
tot_L_agg_nogrant = 100*sum(irf_nogrant.L_agg(1:T_last))/T_last;
tot_Y_corp_nogrant = 100*sum(irf_nogrant.Y_corp(1:T_last))/T_last;
tot_L_corp_nogrant = 100*sum(irf_nogrant.L_corp(1:T_last))/T_last;
tot_L_small_nogrant = 100*sum(irf_nogrant.L_small(1:T_last))/T_last;
tot_output_small_nogrant = 100*sum(irf_nogrant.output_small(1:T_last))/T_last;

% Make table with shocks
FID = fopen(fullfile('tables','tran_shocks.tex'),'w');
fprintf(FID,' \\begin{tabular}{llc} \\hline \n');
fprintf(FID,' Parameter & Description & Value \\\\ \n');
fprintf(FID,' \\hline \n');
fprintf(FID,'$ \\eta_{i}$ & Fraction of impacted small firms & %8.4f \\\\ \n',par.eta_i);
fprintf(FID,'$ \\nu^{n}$ & Productivity shock on impacted firms & %8.4f \\\\ \n',par.v_small);
fprintf(FID,'$ \\nu^{\\lambda}$ & Credit shock on small firms & %8.4f \\\\ \n',par.lambda_shift);
fprintf(FID,'$ \\nu^{M}$ & Shock to the mass of potential entrants & %8.4f \\\\ \n',par.mass_shift);
fprintf(FID,'$ \\nu^{c}$ & Productivity shock the corporate sector & %8.4f \\\\ \n',par.v_corp);
fprintf(FID,'$ \\nu^{d}$ & Preference shock                  & %8.4f \\\\ \n',par.util_shift);
fprintf(FID,'$ \\nu^{l}$ & Labor supply shock                & %8.4f \\\\ \n',par.lsupply_shift);
fprintf(FID,'$ \\rho   $ & Autocorrelation                   & %8.4f \\\\ \n',par.rho_shock);
%fprintf(FID,'$ T       $ & Length of transition              & %d    \\\\ \n',par.T);
fprintf(FID, '\\hline \n \\end{tabular} \n');
fclose(FID);

% Make table with transition moments
FID = fopen(fullfile('tables','transition_moments.tex'),'w');
fprintf(FID,' \\begin{tabular}{lcc} \\hline \\hline \n');
fprintf(FID,' Description & Data & Model      \\\\ \n');
fprintf(FID,' \\hline \n');
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Output, 2020Q2           ",data_mom_trans(1,1),model_mom_trans_baseline(1,1));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Output, 2020Q3           ",data_mom_trans(1,2),model_mom_trans_baseline(1,2));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Consumption, 2020Q2       ",data_mom_trans(2,1),model_mom_trans_baseline(2,1));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Total employment, 2020Q2    ",data_mom_trans(5,1),model_mom_trans_baseline(5,1));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Employment small, 2020Q2     ",data_mom_trans(6,1),model_mom_trans_baseline(6,1));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Exit rate, 2020Q2      ",data_mom_trans(8,1),model_mom_trans_baseline(8,1));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Entry rate, 2020Q2      ",data_mom_trans(10,1),model_mom_trans_baseline(10,1));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Private investment, 2020Q2        ",data_mom_trans(3,1),model_mom_trans_baseline(3,1));
fprintf(FID,"%s & %8.4f & %8.4f   \\\\ \n", "Small firm output, 2020Q2",data_mom_trans(4,1),model_mom_trans_baseline(4,1));
fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);


%% Make table for costs of grants (baseline vs targeted)
clear;close all

%% Set up some useful paths

ResultsDir = fullfile('mat'); %folder where .mat files are stored
% Load results for no grant economy
load(fullfile(ResultsDir,'nogrant.mat'),'agg_tran')
L_small_nogrant = agg_tran.L_small;

% Load results for baseline grant economy
load(fullfile(ResultsDir,'grant_baseline.mat'),'agg_tran','agg')
tot_grant_baseline = agg_tran.tot_grant;
L_small_grant_baseline = agg_tran.L_small;
Y_agg_ss = agg.Y_agg;
L_small_ss = agg.L_small;

% Load results for slim targeted (now called targeted) grant economy
load(fullfile(ResultsDir,'grant_targslim.mat'),'agg_tran')
tot_grant_targslim = agg_tran.tot_grant;
L_small_grant_targslim = agg_tran.L_small;

% Load results for slim targeted large grant economy
load(fullfile(ResultsDir,'grant_targslim_large.mat'),'agg_tran')
tot_grant_targslim_large     = agg_tran.tot_grant;
L_small_grant_targslim_large = agg_tran.L_small;

% Load results for slim targeted small grant economy
load(fullfile(ResultsDir,'grant_targslim_small.mat'),'agg_tran')
tot_grant_targslim_small     = agg_tran.tot_grant;
L_small_grant_targslim_small = agg_tran.L_small;

% employment saved 
t0 = 1;
t1 = 40;

empsaved_grant_baseline = sum(L_small_grant_baseline(t0:t1) - L_small_nogrant(t0:t1))/(t1-t0+1);
empsaved_grant_targslim = sum(L_small_grant_targslim(t0:t1) - L_small_nogrant(t0:t1))/(t1-t0+1);
empsaved_grant_targslim_large = sum(L_small_grant_targslim_large(t0:t1) - L_small_nogrant(t0:t1))/(t1-t0+1);
empsaved_grant_targslim_small = sum(L_small_grant_targslim_small(t0:t1) - L_small_nogrant(t0:t1))/(t1-t0+1);

frac_grant_baseline = tot_grant_baseline/(4*Y_agg_ss);
frac_grant_targslim = tot_grant_targslim/(4*Y_agg_ss);
frac_grant_targslim_large = tot_grant_targslim_large/(4*Y_agg_ss);
frac_grant_targslim_small = tot_grant_targslim_small/(4*Y_agg_ss);

frac_empsaved_baseline = empsaved_grant_baseline/L_small_ss;
frac_empsaved_targslim = empsaved_grant_targslim/L_small_ss;
frac_empsaved_targslim_large = empsaved_grant_targslim_large/L_small_ss;
frac_empsaved_targslim_small = empsaved_grant_targslim_small/L_small_ss;


% Make table with cost of grants (only baseline and targeted
FID = fopen(fullfile('tables','cost_grants.tex'),'w');
fprintf(FID,' \\begin{tabular}{lcc} \\hline \\hline \n');
fprintf(FID,'  & Baseline grant & Targeted grant  \\\\ \n');
fprintf(FID,' \\hline \n');
fprintf(FID,' Cost (Frac. GDP) & %8.4f & %8.4f  \\\\ \n', ...
    frac_grant_baseline,frac_grant_targslim);
fprintf(FID,' Emp. save (Frac. Emp) & %8.4f & %8.4f  \\\\ \n', ...
    frac_empsaved_baseline,frac_empsaved_targslim);
fprintf(FID,' Cost per perc. jobs saved & %8.4f & %8.4f  \\\\ \n', ...
    frac_grant_baseline/(100*frac_empsaved_baseline),frac_grant_targslim/(100*frac_empsaved_targslim));

fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);
% 

% Make table with cost of grants
FID = fopen(fullfile('tables','cost_grants_robust.tex'),'w');
fprintf(FID,' \\begin{tabular}{lcccc} \\hline \\hline \n');
fprintf(FID,'  & Baseline grant & Targeted grant & Targeted grant (large) & Targeted grant (small) \\\\ \n');
fprintf(FID,' \\hline  \n');
fprintf(FID,' Cost (Frac. GDP) & %8.4f & %8.4f & %8.4f & %8.4f \\\\ \n', ...
    frac_grant_baseline,frac_grant_targslim,frac_grant_targslim_large,frac_grant_targslim_small);
fprintf(FID,' Emp. save (Frac. Emp) & %8.4f & %8.4f & %8.4f & %8.4f \\\\ \n', ...
    frac_empsaved_baseline,frac_empsaved_targslim,frac_empsaved_targslim_large,frac_empsaved_targslim_small);
fprintf(FID,' Cost per perc. jobs saved & %8.4f & %8.4f & %8.4f & %8.4f \\\\ \n', ...
    frac_grant_baseline/(100*frac_empsaved_baseline),frac_grant_targslim/(100*frac_empsaved_targslim),...
    frac_grant_targslim_large/(100*frac_empsaved_targslim_large),frac_grant_targslim_small/(100*frac_empsaved_targslim_small));
fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
fclose(FID);

