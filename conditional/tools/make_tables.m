%% Housekeeping

% clear;clc;close all
% format long g
% addpath(genpath('utilities'));

%Some flags, controls etc.
SaveDir                   = 'results\';
par.SaveDir               = SaveDir;
makeCompleteLatexDocument = 1;

%% Make Table VI - Distributional Implications - United States, benchmark

load([SaveDir 'benchmark.mat'])

iTau = iTau_bench; %choose index corresp. to benchmark economy
wealth_for_table = 100*[WealthShares{iTau}.bot40,WealthShares{iTau}.top20,WealthShares{iTau}.top10,WealthShares{iTau}.top1];
earnin_for_table = 100*[EarningsShares{iTau}.bot40,EarningsShares{iTau}.top20,EarningsShares{iTau}.top10,EarningsShares{iTau}.top1];
income_for_table = 100*[IncomeShares{iTau}.bot40,IncomeShares{iTau}.top20,IncomeShares{iTau}.top10,IncomeShares{iTau}.top1];

fid = fopen([par.SaveDir 'tableVI.tex'], 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Distributional Implications - United States} \n');
%fprintf(fid,'\\label{tab:results_us} \n');

fprintf(fid, '\\begin{tabular}{lccccc} \\hline \n');
fprintf(fid, '\\hline \n');
%fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,'                  & Gini & Bottom 40 & Top 20 & Top 10 & Top 1  \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid,'Wealth            &      &           &        &        &       \\\\ \n');
fprintf(fid,'Actual U.S. data  &%8.2f & %8.2f     & %8.2f  & %8.2f  & %8.2f \\\\ \n',0.78,1.4,79.5,66.1,29.5);
fprintf(fid,'Model, benchmark  &%8.2f & %8.2f     & %8.2f  & %8.2f  & %8.2f  \\\\ \n',WealthGini(iTau),wealth_for_table);
fprintf(fid, '\\hline \n');
fprintf(fid,'Earnings          &      &           &        &        &   \\\\ \n');
fprintf(fid,'Actual U.S. data  &%8.2f & %8.2f     & %8.2f  & %8.2f  & %8.2f  \\\\ \n',0.63,2.8,61.4,43.5,14.8);
fprintf(fid,'Model, benchmark  &%8.2f & %8.2f     & %8.2f  & %8.2f  & %8.2f  \\\\ \n',EarningsGini(iTau),earnin_for_table);
fprintf(fid, '\\hline \n');
fprintf(fid,'Total Income      &      &           &        &        &   \\\\ \n');
fprintf(fid,'Actual U.S. data  & %8.2f& %8.2f     & %8.2f  & %8.2f  & %8.2f  \\\\ \n',0.57,8.8,59.9,45.2,18.6);
fprintf(fid,'Model, benchmark  & %8.2f& %8.2f     & %8.2f  & %8.2f  & %8.2f  \\\\ \n',IncomeGini(iTau),income_for_table);
fprintf(fid,'\\hline \n \\end{tabular} \n');
% End tabular environment

fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

%% Make Table VIII, correlations, benchmark

makeCompleteLatexDocument = 1;
fid = fopen([par.SaveDir 'tableVIII.tex'], 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Correlation of Earnings, Income and Wealth} \n');
%fprintf(fid,'\\label{tab:results_us} \n');

fprintf(fid, '\\begin{tabular}{lcc} \\hline \n');
fprintf(fid, '\\hline \n');
%fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,'                       & Data  & Model   \\\\ \n');
fprintf(fid,'                       &       & $\\tau^h=0.36$ (Bench.)  \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid,' corr(\\textit{earnings,income}) & %8.2f & %8.2f  \\\\ \n',0.93,Corr_e_i(iTau_bench));
fprintf(fid,' corr(\\textit{earnings,assets}) & %8.2f & %8.2f  \\\\ \n',0.23,Corr_e_w(iTau_bench));
fprintf(fid,' corr(\\textit{income,assets})   & %8.2f & %8.2f  \\\\ \n',0.32,Corr_i_w(iTau_bench));
fprintf(fid, '\\hline \n');

fprintf(fid,'\\hline \n \\end{tabular} \n');
% End tabular environment

fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

%% Distribution of labor income

income_quintiles = 100*[IncomeShares{iTau_bench}.q1;
                        IncomeShares{iTau_bench}.q2;
                        IncomeShares{iTau_bench}.q3;
                        IncomeShares{iTau_bench}.q4;
                        IncomeShares{iTau_bench}.q5];
                    
income_top = 100*[IncomeShares{iTau_bench}.top90_95;
                  IncomeShares{iTau_bench}.top5;
                  IncomeShares{iTau_bench}.top1];

fid = fopen([par.SaveDir 'table_income.tex'], 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Shares of (labor) income - model and data} \n');
%fprintf(fid,'\\label{tab:results_us} \n');

fprintf(fid, '\\begin{tabular}{lcc} \\hline \n');
fprintf(fid, '\\hline \n');
%fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,' Percentiles of income & Data & Model \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid,' \\textit{Quantile} &  &  \\\\ \n');
fprintf(fid,'1st (bottom 20%%) & %8.2f & %8.2f     \\\\ \n',2.1,   income_quintiles(1));
fprintf(fid,'2nd (20-40\\%%)     & %8.2f & %8.2f     \\\\ \n',6.7, income_quintiles(2));
fprintf(fid,'3rd (40-60\\%%)     & %8.2f & %8.2f     \\\\ \n',12.3,income_quintiles(3));
fprintf(fid,'4th (60-80\\%%)     & %8.2f & %8.2f     \\\\ \n',21.3,income_quintiles(4));
fprintf(fid,'5th (80-100\\%%)    & %8.2f & %8.2f     \\\\ \n',57.6,income_quintiles(5));
fprintf(fid,' \\textit{Top}  &  &  \\\\ \n');
fprintf(fid,'90-95\\%%           & %8.2f & %8.2f     \\\\ \n',11.7,income_top(1));
fprintf(fid,'5\\%%               & %8.2f & %8.2f     \\\\ \n',29.1,income_top(2));
fprintf(fid,'1\\%%               & %8.2f & %8.2f     \\\\ \n',14.3,income_top(3));
fprintf(fid, '\\hline \n');
fprintf(fid,'Gini coefficient & %8.2f & %8.2f    \\\\ \n',0.55,IncomeGini(iTau_bench));

fprintf(fid,'\\hline \n \\end{tabular} \n');
% End tabular environment

fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

%% Distribution of taxes paid, by percentiles of income
%Consider benchmark economy only
taxes_quintiles = 100*taxes_quintiles_mat(iTau_bench,:);

fid = fopen([par.SaveDir 'table_taxes.tex'], 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Shares of tax payments - model and data} \n');
%fprintf(fid,'\\label{tab:results_us} \n');

fprintf(fid, '\\begin{tabular}{lcc} \\hline \n');
fprintf(fid, '\\hline \n');
%fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,' Percentiles of income & Data & Model \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid,' \\textit{Quantile} &  &  \\\\ \n');
fprintf(fid,'1st (bottom 20\\%%) & %8.2f & %8.2f     \\\\ \n',0.3,taxes_quintiles(1));
fprintf(fid,'2nd (20-40\\%%)     & %8.2f & %8.2f     \\\\ \n',2.2,taxes_quintiles(2));
fprintf(fid,'3rd (40-60\\%%)     & %8.2f & %8.2f     \\\\ \n',6.9,taxes_quintiles(3));
fprintf(fid,'4th (60-80\\%%)     & %8.2f & %8.2f     \\\\ \n',15.9,taxes_quintiles(4));
fprintf(fid,'5th (80-100\\%%)    & %8.2f & %8.2f     \\\\ \n',74.6,taxes_quintiles(5));

fprintf(fid,'\\hline \n \\end{tabular} \n');
% End tabular environment

fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

stop

%% Make Table IX, results, benchmark + experiments

load([SaveDir 'experiments.mat'])

makeCompleteLatexDocument = 1;

benchmark_vec = [tau_h_vec(iTau_bench),Welfare_vec(iTau_bench),r(iTau_bench)*100,K_Y_vec(iTau_bench),H_vec(iTau_bench),Y(iTau_bench),...
    C_Y_vec(iTau_bench),hours_vec(iTau_bench),B_Y_vec(iTau_bench),T_Y_vec(iTau_bench),T_vec(iTau_bench) ]; %row vector

%Experiments are computed for many tax rates.. wanna select only a subset
chosen = (0.25:0.05:0.65)';
indVec = zeros(length(chosen),1);
for i=1:length(chosen)
    indVec(i,1) = find(tau_h_vec==chosen(i));
end

experiments_vec = zeros(length(chosen),11);
for ii=1:length(chosen)
    iTau = indVec(ii);
    experiments_vec(ii,1:11) = [tau_h_vec(iTau),Welfare_vec(iTau),r(iTau)*100,K_Y_vec(iTau),H_vec(iTau),Y(iTau),...
    C_Y_vec(iTau),hours_vec(iTau),B_Y_vec(iTau),T_Y_vec(iTau),T_vec(iTau) ]; %row vector
end

fid = fopen([par.SaveDir 'tableIX.tex'], 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Results for Different Tax Rates - United States} \n');
fprintf(fid,'\\label{tab:results_us} \n');

fprintf(fid, '\\begin{tabular}{ccccccccccc} \\hline \n');
fprintf(fid, '\\hline \n');
%fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,'$\\tau^h$ & w.g. & $r$ & $K/Y$ & $H$ & $Y$ & $C/Y$ & $h$ & $B/Y$ & $T/Y$ & $T$ \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid,'Benchmark &  &  &  &  &  &  &  &  &  &  \\\\ \n');
fprintf(fid,'%8.3f & %8.2f & %8.2f & %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',benchmark_vec);
fprintf(fid,'Experiments &  &  &  &  &  &  &  &  &  &  \\\\ \n');
for ii=1:length(chosen)
    fprintf(fid,'%8.2f & %8.2f & %8.2f & %8.2f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',experiments_vec(ii,:));
end
fprintf(fid,'\\hline \n \\end{tabular} \n');
% End tabular environment

fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);






%fprintf(fid,'Avg. hours worked (\\%%) & %8.3f & %8.3f \\\\ \n',33 ,(ave_n_work/maxLwork)*100 );

%% Make and save some plots

% [junk,ind_max] = max(WelfareVec);
% tau_n_max      = tau_h_vec(ind_max);
% 
% figure(1)
% plot(tau_h_vec,WelfareVec,'linewidth',2)
% hold on
% plot(tau_n_max,WelfareVec(ind_max),'o','linewidth',4)
% grid on
% xlabel('Tax Rate on Labor Income, \tau^{n}','fontsize',14)
% ylabel('Aggregate Welfare','fontsize',14)
% title('Welfare-maximizing Labor Tax','fontsize',14)
% print([SaveDir 'welfare_tax_n'],'-dpng')
% 
% %% Write results of the calibrations in latex tables
% 
% %% Table 4
% % Each cell entry to each calibration
% calib = cell(length(tau_h_vec),1);
% for i=1:length(tau_h_vec)
%     calib{i} = [WelfareVec(i), Corr_h_z(i), HoursCV(i), H(i), K(i)/Y(i), w(i)*L(i)/Y(i), I(i)/Y(i), TaxesToY(i) ]';
% end
% 
% varNames = {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10'}; %cell array of strings
% rowNames = {'Welfare','corr(h,z)','cv(h)','H','K/Y','wL/Y','I/Y','Taxes/Y'}; %cell array of strings
% T4 = table(calib{1},calib{2},calib{3},calib{4},calib{5},calib{6},calib{7},...
%     calib{8},calib{9},calib{10},'VariableNames',varNames,'RowNames',rowNames);
% % Apply func (here, round) to all entries of table T1
% T4 = varfun(@(x) round(x, 3), T4);
% T4.Properties.VariableNames = varNames;
% T4.Properties.RowNames      = rowNames;
% disp(T4)
% 
% scaleBox = 1;
% makeFullTexDoc = 1;
% makeTableFloat = 1;
% table2latexV2(T4, [SaveDir 'table4'],scaleBox,makeFullTexDoc,makeTableFloat);

% %% Table 5
% % Each cell entry to each calibration
% calib = cell(4,1);
% for i=1:4
%     calib{i} = [Y(i) C(i) K(i) L(i) H(i) L(i)/H(i) Y(i)/H(i) Corr_h_z(i) Corr_k_z(i)]';
% end
% 
% varNames = {'C1','C2','C3','C4'}; %cell array of strings
% rowNames = {'Y','C','K','L','H','L/H','Y/H','corr(h,z)','corr(k,z)'}; %cell array of strings
% 
% T1 = table(calib{1},calib{2},calib{3},calib{4},'VariableNames',varNames,'RowNames',rowNames);
% % Apply func (here, round) to all entries of table T1
% T1 = varfun(@(x) round(x, 3), T1);
% T1.Properties.VariableNames = varNames;
% T1.Properties.RowNames      = rowNames;
% disp(T1)
% 
% %Unfortunately, writetable does not accept .tex but only txt and excel
% writetable(T1,[SaveDir 'table5.txt']);
% writetable(T1,[SaveDir 'table5.xlsx'],'WriteRowNames',true);
% 
% %This writes table directly into a latex file. Unfortunately it's not a
% %stand-alone latex document. It automatically writes row names and variable names. 
% full_latex = 1;
% table2latexV2(T1, [SaveDir 'table5'],full_latex);
% 
% %% Table 6
% 
% varNames = {'cv','Gini','q1','q2','q3','q4','q5'}; %cell array of strings
% rowNames = {'Earnings C1';
%             'Earnings C2';
%             'Earnings C3';
%             'Earnings C4';
%             'Income C1';
%             'Income C2';
%             'Income C3';
%             'Income C4';
%             'Wealth C1';
%             'Wealth C2';
%             'Wealth C3';
%             'Wealth C4'};
% 
% T6 = table([EarningsCV;IncomeCV;WealthCV] ,[EarningsGini;IncomeGini;WealthGini],...
%     [EarningsQuintileShares(:,1);IncomeQuintileShares(:,1);WealthQuintileShares(:,1)],...
%     [EarningsQuintileShares(:,2);IncomeQuintileShares(:,2);WealthQuintileShares(:,2)],...
%     [EarningsQuintileShares(:,3);IncomeQuintileShares(:,3);WealthQuintileShares(:,3)],...
%     [EarningsQuintileShares(:,4);IncomeQuintileShares(:,4);WealthQuintileShares(:,4)],...
%     [EarningsQuintileShares(:,5);IncomeQuintileShares(:,5);WealthQuintileShares(:,5)],...
%     'VariableNames',varNames,'RowNames',rowNames);
% 
% T6 = varfun(@(x) round(x, 2), T6);
% T6.Properties.VariableNames = varNames;
% T6.Properties.RowNames      = rowNames;
% 
% disp(T6);
% 
% full_latex = 1;
% table2latexV2(T6, [SaveDir 'table6'],full_latex);