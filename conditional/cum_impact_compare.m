%% Compare short-run vs long-run cumulative effect of agg shock under different policy environments
% The script loads the following mat files:
% - nogrant.mat
% - grant_baseline.mat
% - grant_targslim.mat

clear;clc;close all

matNames = {'nogrant','grant_baseline'};
for ii = 1:numel(matNames)
    if ~isfile(fullfile('mat',[matNames{ii},'.mat']))
        warning('MAT file ''%s'' is missing \n',matNames{ii})
        pause
    end
end

disp('Cumulative impacts of the pandemic')

%% Set up some useful paths

ResultsDir = fullfile('mat'); %folder where .mat files are stored
FormatFig = '-depsc';  % Specify '-dpng' or '-depsc'
% Specify here where you want to save the figures

%% Make IRFs for no grant economy
disp('NO GRANT ECONOMY:')
% Load results for no grant economy
if isfile(fullfile('mat','nogrant.mat'))
    load(fullfile('mat','nogrant.mat'),'irf')
else
    error('File nogrant.mat does not exist!')
end

irf_nogrant=irf;
clear irf

%% Make IRFs for baseline grant economy
disp('BASELINE GRANT ECONOMY:')
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'irf')
else
    error('File grant_baseline.mat does not exist!')
end

% Compute IRFs in percentage deviation from steady-state
irf_grant_baseline = irf;
clear irf

%% Compare cumulative impact: Short run, medium run, long run

% Number of quarters in the short-run, mid-run, and long-run
T_sr = 1; % 2 quarters % we changed T_sr from 2 quaters to 1q
T_mr = 8; % 3 years
T_lr = 40; % 10 years


varNames = {'C_agg','K_agg','K_small','K_corp','mass_small',...
    'liq','exit_rate','exit','entry_rate','Y_agg','Y_corp','output_small','L_agg','L_small','L_corp',...
    'KL_ratio','InvK_corp','InvK','q','w'};

varNamesLabels = {'Consumption';'Agg. capital';'Capital, small firms';'Capital, corp.';'Mass of small firms';
    'Liquidation value';'Exit rate';'Measure of exiting firms'; 'Entry rate'; 'Agg. output'; 'Output, corp.';'Output, small firms';
    'Agg. emp.';'Emp., small firms';'Emp., corp.';'K-L ratio';'Investment, corp.';'Agg. investment'
    'Price of financial asset';'Wage'};


% Create a stacked bar chart using the bar function
for ii = 1:numel(varNames)
    varii = varNames{ii};
    sr_totimp_nogrant.(varii) = 100*sum(irf_nogrant.(varii)(1:T_sr))/T_lr;
    sr_totimp_grant_baseline.(varii) = 100*sum(irf_grant_baseline.(varii)(1:T_sr))/T_lr;
     
    mr_totimp_nogrant.(varii) = 100*sum(irf_nogrant.(varii)(T_sr+1:T_mr))/T_lr;
    mr_totimp_grant_baseline.(varii) = 100*sum(irf_grant_baseline.(varii)(T_sr+1:T_mr))/T_lr;
     
    lr_totimp_nogrant.(varii) = 100*sum(irf_nogrant.(varii)(T_mr+1:T_lr))/T_lr;
    lr_totimp_grant_baseline.(varii) = 100*sum(irf_grant_baseline.(varii)(T_mr+1:T_lr))/T_lr;
     
    net_totimp_nogrant.(varii) = sr_totimp_nogrant.(varii) + mr_totimp_nogrant.(varii) +lr_totimp_nogrant.(varii) ;
    net_totimp_grant_baseline.(varii) = sr_totimp_grant_baseline.(varii) + mr_totimp_grant_baseline.(varii) +lr_totimp_grant_baseline.(varii) ;
 

end

%% Make Stacked Bar Charts
disp('Make Stacked Bar Charts...')
FigDir = fullfile('figures','grant_vs_nogrant'); %folder where output figures are stored
myfontsize = 20;

% Only baseline grant vs no grant
for ii = 1:numel(varNames)
   
    varii = varNames{ii};
    varlabel =  varNamesLabels{ii};
    figname = ['bar_totimp_',varii];
    sr_totimp.(varii) = [sr_totimp_nogrant.(varii) sr_totimp_grant_baseline.(varii) ]';
    mr_totimp.(varii) = [mr_totimp_nogrant.(varii)   mr_totimp_grant_baseline.(varii)]';
    lr_totimp.(varii) = [lr_totimp_nogrant.(varii)  lr_totimp_grant_baseline.(varii) ]';
    net_totimp.(varii) = [net_totimp_nogrant.(varii)  net_totimp_grant_baseline.(varii) ]';
    
    figure('Position',[1000 918 560 420])
    b=bar(1:length(sr_totimp.(varii) ), [sr_totimp.(varii)  mr_totimp.(varii)  lr_totimp.(varii) ], 0.5, 'stack');
    hold on
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Conditional grant'},'fontsize',myfontsize)
    xtickangle(30)
    hL=plot(1:length(net_totimp.(varii)),net_totimp.(varii),'sg','MarkerSize',15,'MarkerFaceColor','g');
    hold off
    %title(varlabel,'fontsize',myfontsize)
    ylabel('% Impact','fontsize',myfontsize)
    ax = gca;
    ax.YAxis.Exponent = 0;
    legend('Q1', 'Q2-Q8', 'Q9-Q40','Average','location','best','fontsize',myfontsize-2)
    legend('boxoff')
    print(fullfile(FigDir,figname),FormatFig);
end




