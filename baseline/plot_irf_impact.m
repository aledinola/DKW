function [] = plot_irf_impact(irf_impact,SaveDir,do_save,FormatFig)

% INPUTS:
%   "irf_impact"  Structure with IRF by impact status
%   "SaveDir"     Folder where to save the plots
%   "do_save"     Flag 0/1 to save plots or not 
%   "FormatFig"   '-dpng' or '-depsc' 

% Check inputs
if ~isstruct(irf_impact)
    error('Input irf_impact in plot_irf_impact must be a structure')
end

%% Options for plots
lastp      = 20;%32;
myfontsize = 20;
mylw       = 4;

%% Compute IRF differences
varNames = {'K_small_owned','mass_small','exit_rate','output_small','L_small'};

varNamesLabels = {'Capital, small firms';'Mass of small firms';
    'Exit rate'; 'Output, small firms';'Emp., small firms'};


%% Plot irf
for ii=1:numel(varNames)
    name = varNames{ii};
    label = varNamesLabels{ii};

    % Plot 
    figure
    plot(irf_impact.(name)(1:lastp,1)*100,'-',"linewidth",mylw)
    hold on
    plot(irf_impact.(name)(1:lastp,2)*100,'-',"linewidth",mylw)
    hold on
    yline(0,'--');
    legend('Impacted','Not impacted','FontSize',myfontsize,'Location','Best')
    xlabel("Time in transition, t",'FontSize',myfontsize)
    title(label,'FontSize',myfontsize)
    figname = ['irf_impact_',name];
    if do_save==1; print(fullfile(SaveDir,figname),FormatFig); end
    

end



end % end function plot_irf_impact

