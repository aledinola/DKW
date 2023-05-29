function [] = plot_ppchange_compare(irf1,irf2,name1,name2,lastp,SaveDir,do_save,FormatFig)

% DESCRIPTION
%   plot_ppchange_compare plots pctage point change for two different scenarios 
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
varNames = {'exit_rate','entry_rate'};

varNamesLabels = {'Exit rate'; 'Entry rate'};
   
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
    ylabel("Percentage point change",'FontSize',myfontsize)
    grid on
    %title(label,'FontSize',myfontsize)
    figname = ['ppchange_',name];
    if do_save==1; print(fullfile(SaveDir,figname),FormatFig); end

end



end % end function plot_ppchange_compare

