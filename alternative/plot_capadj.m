function [] = plot_capadj(capadj_rate_ss,capadj_rate_baseline,capadj_rate_nogrant,capadj_rate_targslim,...
    sr_cum_capadj_nogrant,sr_cum_capadj_grant_baseline,sr_cum_capadj_grant_targslim,...
    sum_sr_cum_capadj_nogrant,sum_sr_cum_capadj_grant_baseline,sum_sr_cum_capadj_grant_targslim,...
    mr_cum_capadj_nogrant,mr_cum_capadj_grant_baseline,mr_cum_capadj_grant_targslim,...
    sum_mr_cum_capadj_nogrant,sum_mr_cum_capadj_grant_baseline,sum_mr_cum_capadj_grant_targslim,...
    lr_cum_capadj_nogrant,lr_cum_capadj_grant_baseline,lr_cum_capadj_grant_targslim,...
    sum_lr_cum_capadj_nogrant,sum_lr_cum_capadj_grant_baseline,sum_lr_cum_capadj_grant_targslim,...
    ave_cum_capadj_nogrant,ave_cum_capadj_grant_baseline,ave_cum_capadj_grant_targslim,...
    sum_ave_cum_capadj_nogrant,sum_ave_cum_capadj_grant_baseline,sum_ave_cum_capadj_grant_targslim,...
    T_last,LW,FS,SaveDir,fileroot,FormatFig,do_save)
% plot_capadj plots capital adjustment in the transition under the baseline
% and counterfactual scenarios.

%%
figure
plot(0:T_last-1,100*capadj_rate_baseline(1:T_last,1),'--','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_nogrant(1:T_last,1),'-','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_targslim(1:T_last,1),':','LineWidth',LW)
yline(100*capadj_rate_ss(1),'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% of steady state capital in small firms','FontSize',FS)
legend('Baseline','Laissez-Faire','Targeted','Steady State','FontSize',FS,...
    'location','best')
title('Upward capital adjustment by active firms')
set(gca,'FontSize',FS)
filename = ['capadj_diff_1',fileroot];
if do_save==1; print(fullfile(SaveDir,filename),FormatFig); end

figure
plot(0:T_last-1,100*capadj_rate_baseline(1:T_last,2),'--','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_nogrant(1:T_last,2),'-','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_targslim(1:T_last,2),':','LineWidth',LW)
yline(100*capadj_rate_ss(2),'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% of steady state capital in small firms','FontSize',FS)
legend('Baseline','Laissez-Faire','Targeted','Steady State','FontSize',FS,'location','best')
title('Downward capital adjustment by active firms')
set(gca,'FontSize',FS)
filename = ['capadj_diff_2',fileroot];
if do_save==1; print(fullfile(SaveDir,filename),FormatFig); end

figure
plot(0:T_last-1,100*capadj_rate_baseline(1:T_last,3),'--','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_nogrant(1:T_last,3),'-','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_targslim(1:T_last,3),':','LineWidth',LW)
yline(100*capadj_rate_ss(3),'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% of steady state capital in small firms','FontSize',FS)
legend('Baseline','Laissez-Faire','Targeted','Steady State','FontSize',FS,'location','best')
title('Capital bought by entering firms')
set(gca,'FontSize',FS)
filename = ['capadj_diff_3',fileroot];
if do_save==1; print(fullfile(SaveDir,filename),FormatFig); end

figure
plot(0:T_last-1,100*capadj_rate_baseline(1:T_last,4),'--','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_nogrant(1:T_last,4),'-','LineWidth',LW)
hold on
plot(0:T_last-1,100*capadj_rate_targslim(1:T_last,4),':','LineWidth',LW)
yline(100*capadj_rate_ss(4),'--')
axis tight
xlabel('Time in transition, t','FontSize',FS)
ylabel('% of steady state capital in small firms','FontSize',FS)
legend('Baseline','Laissez-Faire','Targeted','Steady State','FontSize',FS,'location','best')
title('Capital sold by exiting firms')
set(gca,'FontSize',FS)
filename = ['capadj_diff_4',fileroot];
if do_save==1; print(fullfile(SaveDir,filename),FormatFig); end


%% Cumulative capital adjstment in the short run: t = 1 to T_sr
figure('Position',[1000 918 560 420])
b=bar(1:3, [sr_cum_capadj_nogrant, sr_cum_capadj_grant_baseline, sr_cum_capadj_grant_targslim]', 0.5, 'stack');
hold on
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Baseline grant' ...
    'Targeted grant' },'fontsize',FS)
xtickangle(30)
hL=plot(1:3,[sum_sr_cum_capadj_nogrant sum_sr_cum_capadj_grant_baseline ...
    sum_sr_cum_capadj_grant_targslim],'sg','MarkerSize',15,'MarkerFaceColor','g');
hold off
%title(varlabel,'fontsize',myfontsize)
ylabel({'Cumulative excess investment', '(% of steady state capital)'},'fontsize',FS)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = FS;
legend('Invest', 'Disinvest','Entry', 'Exit','All adj.','location','northoutside',...
    'fontsize',FS-1,'NumColumns',5)
filename = ['cum_capadj_decomp_sr',fileroot];
print(fullfile(SaveDir,filename),FormatFig);

% Cumulative capital adjstment in the medium run: t = T_sr+1 to T_mr
figure('Position',[1000 918 560 420])
b=bar(1:3, [mr_cum_capadj_nogrant, mr_cum_capadj_grant_baseline, mr_cum_capadj_grant_targslim]', 0.5, 'stack');
hold on
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Baseline grant' ...
    'Targeted grant' },'fontsize',FS)
xtickangle(30)
hL=plot(1:3,[sum_mr_cum_capadj_nogrant sum_mr_cum_capadj_grant_baseline ...
    sum_mr_cum_capadj_grant_targslim],'sg','MarkerSize',15,'MarkerFaceColor','g');
hold off
%title(varlabel,'fontsize',myfontsize)
ylabel({'Cumulative excess investment', '(% of steady state capital)'},'fontsize',FS)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = FS;
legend('Invest', 'Disinvest','Entry', 'Exit','All adj.','location','northoutside',...
    'fontsize',FS-1,'NumColumns',5)
filename = ['cum_capadj_decomp_mr',fileroot];
print(fullfile(SaveDir,filename),FormatFig);

% Cumulative capital adjstment in the long run: t = T_mr+1 to T_lr
figure('Position',[1000 918 560 420])
b=bar(1:3, [lr_cum_capadj_nogrant, lr_cum_capadj_grant_baseline, lr_cum_capadj_grant_targslim]', 0.5, 'stack');
hold on
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Laissez-Faire' 'Baseline grant' ...
    'Targeted grant' },'fontsize',FS)
xtickangle(30)
hL=plot(1:3,[sum_lr_cum_capadj_nogrant sum_lr_cum_capadj_grant_baseline ...
    sum_lr_cum_capadj_grant_targslim],'sg','MarkerSize',15,'MarkerFaceColor','g');
hold off
%title(varlabel,'fontsize',myfontsize)
ylabel({'Cumulative excess investment', '(% of steady state capital in small firms)'},'fontsize',FS)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = FS;
legend('Invest', 'Disinvest','Entry', 'Exit','All adj.','location','northoutside',...
    'fontsize',FS-1,'NumColumns',5)
filename = ['cum_capadj_decomp_lr',fileroot];
print(fullfile(SaveDir,filename),FormatFig);

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
filename = ['cum_capadj_decomp_all',fileroot];
print(fullfile(SaveDir,filename),FormatFig);

end