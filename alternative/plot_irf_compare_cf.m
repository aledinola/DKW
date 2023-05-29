function [] = plot_irf_compare_cf(irf1,irf2,irf_nogrant,name1,name2,name1_short,name2_short,data_mom_trans,SaveDir,do_save,FormatFig)

% INPUTS:
%   "irf1"        Structure with IRF counterfactual 1 
%   "irf2"        Structure with IRF counterfactual 1 
%   "irf_nogrant" Structure with IRF for no grant economy
%   "name1"       String for legend referring to irf1
%   "name2"       String for legend referring to irf2
%   "data_mom_trans"
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
if ~isstruct(irf_nogrant)
    error('Input irf_nogrant in plot_irf_compare must be a structure')
end

lastp = 32;
myfontsize = 20;
mylw = 4;

varNames = {'C_agg','K_agg','K_all','K_corp','K_small','K_small_owned','mass_small','exit',...
    'entry','Y_agg','Y_corp','output_small','L_agg','L_small','L_corp',...
    'KL_ratio','InvK','InvK_all','q','w','rental'};

% Compute IRF|1 minus IRF|no grant
% Compute IRF|2 minus IRF|no grant
dirf1 = struct();
dirf2 = struct();
for ii=1:numel(varNames)
    name = varNames{ii};
    dirf1.(name) = irf1.(name)-irf_nogrant.(name);
    dirf2.(name) = irf2.(name)-irf_nogrant.(name);
end

%% Plot IRF differences

figure
plot(dirf1.C_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.C_agg(1:lastp)*100,'-.',"linewidth",mylw)
hold on
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Consumption','FontSize',myfontsize)
plotname = ['C_agg_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.K_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.K_agg(1:lastp)*100,'-.',"linewidth",mylw)
hold on
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Household Capital','FontSize',myfontsize)
plotname = ['K_agg_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.K_all(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.K_all(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total capital in the economy','FontSize',myfontsize)
plotname = ['K_all_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end


figure
plot(dirf1.K_corp(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.K_corp(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Capital in corporate sector','FontSize',myfontsize)
plotname = ['K_corp_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end


figure
plot(dirf1.K_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.K_small(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('(Rented) capital in small firms','FontSize',myfontsize)
plotname = ['K_small_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.K_small_owned(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.K_small_owned(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Owned capital in small firms','FontSize',myfontsize)
plotname = ['K_small_owned_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.mass_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.mass_small(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Mass of small firms','FontSize',myfontsize)
plotname = ['mass_small_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.exit(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.exit(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Exit','FontSize',myfontsize)
plotname = ['exit_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.entry(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.entry(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Entry','FontSize',myfontsize)
plotname = ['entry_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.Y_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.Y_agg(1:lastp)*100,'-.',"linewidth",mylw)
hold on
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total output (small firms + corporate)','FontSize',myfontsize)
plotname = ['Y_agg_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.Y_corp(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.Y_corp(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Output corporate sector','FontSize',myfontsize)
plotname = ['Y_corp_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.output_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.output_small(1:lastp)*100,'-.',"linewidth",mylw)
hold on
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Output small firms','FontSize',myfontsize)
plotname = ['output_small_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end


figure
plot(dirf1.L_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.L_agg(1:lastp)*100,'-.',"linewidth",mylw)
hold on
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total Employment','FontSize',myfontsize)
plotname = ['L_agg_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end


figure
plot(dirf1.L_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.L_small(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Employment in small firms','FontSize',myfontsize)
plotname = ['L_small_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end


figure
plot(dirf1.L_corp(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.L_corp(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Employment in the corporate sector','FontSize',myfontsize)
plotname = ['L_corp_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.KL_ratio(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.KL_ratio(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Capital-to-Labor ratio','FontSize',myfontsize)
plotname = ['KL_ratio_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.InvK(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.InvK(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Household Investment','FontSize',myfontsize)
plotname = ['InvK_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.InvK_all(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.InvK_all(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total Investment in the economy','FontSize',myfontsize)
plotname = ['InvK_all_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.q(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.q(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Financial discount factor','FontSize',myfontsize)
plotname = ['q_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.w(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.w(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Wage','FontSize',myfontsize)
plotname = ['w_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(dirf1.rental(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(dirf2.rental(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Rental rate','FontSize',myfontsize)
plotname = ['rental_dirf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

%% IRF plots
figure
plot(irf1.C_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf2.C_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.C_agg(1:lastp)*100,'-.',"linewidth",mylw)
hold on
%plot(data_mom_trans(2,1:3),'--',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Consumption, % change','FontSize',myfontsize)
plotname = ['C_agg_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(irf1.K_agg(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.K_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.K_agg(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Household Capital, % change','FontSize',myfontsize)
plotname = ['K_agg_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.K_all(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.K_all(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.K_all(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total capital in the economy, % change','FontSize',myfontsize)
plotname = ['K_all_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(irf1.K_corp(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.K_corp(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.K_corp(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Capital in corporate sector, % change','FontSize',myfontsize)
plotname = ['K_corp_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.K_small(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.K_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.K_small(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Rented capital in small firms, % change','FontSize',myfontsize)
plotname = ['K_small_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.K_small_owned(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.K_small_owned(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.K_small_owned(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Owned capital in small firms, % change','FontSize',myfontsize)
plotname = ['K_small_owned_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(irf1.mass_small(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.mass_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.mass_small(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Mass of small firms, % change','FontSize',myfontsize)
plotname = ['mass_small_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.exit(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.exit(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.exit(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Exit, % change','FontSize',myfontsize)
plotname = ['exit_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.entry(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.entry(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.entry(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Entry, % change','FontSize',myfontsize)
plotname = ['entry_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.Y_agg(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.Y_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.Y_agg(1:lastp)*100,'-.',"linewidth",mylw)
hold on
%plot(data_mom_trans(1,1:3),'--',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total output (small firms + corporate), % change','FontSize',myfontsize)
plotname = ['Y_agg_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.Y_corp(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.Y_corp(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.Y_corp(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Output corporate sector, % change','FontSize',myfontsize)
plotname = ['Y_corp_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.output_small(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.output_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.output_small(1:lastp)*100,'-.',"linewidth",mylw)
hold on
%plot(data_mom_trans(4,1:3),'--',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Output small firms, % change','FontSize',myfontsize)
plotname = ['output_small_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end




figure
plot(irf1.L_agg(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.L_agg(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.L_agg(1:lastp)*100,'-.',"linewidth",mylw)
hold on
%plot(data_mom_trans(5,1:4),'--',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total Employment, % change','FontSize',myfontsize)
plotname = ['L_agg_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end




figure
plot(irf1.L_small(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.L_small(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.L_small(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Employment in small firms, % change','FontSize',myfontsize)
plotname = ['L_small_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end




figure
plot(irf1.L_corp(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.L_corp(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.L_corp(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Employment in the corporate sector, % change','FontSize',myfontsize)
plotname = ['L_corp_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end


figure
plot(irf1.KL_ratio(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.KL_ratio(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.KL_ratio(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Capital-to-Labor ratio, % change','FontSize',myfontsize)
plotname = ['KL_ratio_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.InvK(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.InvK(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.InvK(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Household Investment, % change','FontSize',myfontsize)
plotname = ['InvK_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.InvK_all(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.InvK_all(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.InvK_all(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Total Investment in the economy, % change','FontSize',myfontsize)
plotname = ['InvK_all_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(irf1.q(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.q(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.q(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Financial discount factor, % change','FontSize',myfontsize)
plotname = ['q_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end

figure
plot(irf1.w(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.w(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.w(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Wage, % change','FontSize',myfontsize)
plotname = ['w_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



figure
plot(irf1.rental(1:lastp)*100,':',"linewidth",mylw)
hold on
plot(irf2.rental(1:lastp)*100,'-',"linewidth",mylw)
hold on
plot(irf_nogrant.rental(1:lastp)*100,'-.',"linewidth",mylw)
yline(0,'--');
legend(name1,name2,'Laissez-faire','FontSize',myfontsize,'Location','Best')
xlabel("Time in transition, t",'FontSize',myfontsize)
ylabel('Rental rate, % change','FontSize',myfontsize)
plotname = ['rental_irf_',name1_short,'_',name2_short];
if do_save==1; print(fullfile(SaveDir,plotname),FormatFig); end



end % end function plot_irf_compare

