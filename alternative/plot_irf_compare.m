function [] = plot_irf_compare(irf1,irf2,name1,name2,lastp,SaveDir,do_save,FormatFig)

% DESCRIPTION
%   plot_irf_compare plots IRFs for two different scenarios 
% INPUTS:
%   "irf1"        Structure with IRF scenario 1 
%   "irf2"        Structure with IRF scenario 2 
%   "name1"       String for legend referring to irf1
%   "name2"       String for legend referring to irf2
%   "SaveDir"     Folder where to save the plots
%   "do_save"     Flag 0/1 to save plots or not 
%   "FormatFig"   '-dpng' or '-depsc' 

% Check inputs
if ~ischar(name1)
    error('Input name1 in plot_irf_compare must be a string')
end
if ~ischar(name2)
    error('Input name2 in plot_irf_compare must be a string')
end
if ~isstruct(irf1)
    error('Input irf1 in plot_irf_compare must be a structure')
end
if ~isstruct(irf2)
    error('Input irf2 in plot_irf_compare must be a structure')
end

%% Options for plots
%lastp      = length(irf1.C_agg);%32;
myfontsize = 20;
mylw       = 4;

%% Compute IRF differences
varNames = {'C_agg','K_agg','K_small','K_corp','mass_small',...
    'liq','exit_rate','exit','entry_rate','Y_agg','Y_corp','output_small','L_agg','L_small','L_corp',...
    'KL_ratio','InvK_corp','InvK','q','w'};

varNamesLabels = {'Consumption';'Agg. capital';'Capital, small firms';'Capital, corp.';'Mass of small firms';
    'Liquidation value';'Exit rate';'Measure of exiting firms'; 'Entry rate'; 'Agg. output'; 'Output, corp.';'Output, small firms';
    'Agg. emp.';'Emp., small firms';'Emp., corp.';'K-L ratio';'Investment, corp.';'Agg. investment'
    'Price of financial asset';'Wage'};
   
if numel(varNames)~=numel(varNamesLabels)
    error('varNames and varNamesLabels must have the same no. of elements')
end

%% Plot dirf and irf

for ii=1:numel(varNames)
    name = varNames{ii};
    label = varNamesLabels{ii};

    % Plot irf
    figure
    plot(0:lastp-1,irf1.(name)(1:lastp)*100,'-',"linewidth",mylw)
    hold on
    plot(0:lastp-1,irf2.(name)(1:lastp)*100,'--o',"linewidth",mylw-1)
    yline(0,'--');
    legend(name1,name2,'FontSize',myfontsize,'Location','Best')
    ax = gca;
    ax.FontSize = myfontsize;
    xlabel("Time in transition, t",'FontSize',myfontsize)
    ylabel("% change",'FontSize',myfontsize)
    grid on
    %title(label,'FontSize',myfontsize)
    figname = ['irf_',name];
    if do_save==1; print(fullfile(SaveDir,figname),FormatFig); end

end



end % end function plot_irf_compare

