clear
clc
close all

load data

eta_vec = [0.5, 1]; % for comparison 

%% results for eta = 0.5
agg_tran = results(1).agg_tran;
path     = results(1).path;

figure
plot(agg_tran.b_grid(1,:,1),'-o',"linewidth",2)
hold on
plot(agg_tran.b_grid(1,:,2),'-o',"linewidth",2)
yline(b_grid(1),'--',{'Steady-state'})
xlabel("Time in transition, t")
title("b tilde, first point of dynamic b grid")
legend("Impacted","Unimpacted")
hold off
if do_save==1; print(fullfile(FigDir,'b_grid'),'-dpng'); end

figure
plot(agg_tran.K_agg,"linewidth",2)
yline(agg.K_agg,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Aggregate capital K")
if do_save==1; print(fullfile(FigDir,'K'),'-dpng'); end

figure
plot(path.C,"linewidth",2)
yline(agg.C_agg,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Consumption")
if do_save==1; print(fullfile(FigDir,'C'),'-dpng'); end

figure
plot(path.KL_ratio,"linewidth",2)
yline(agg.K_corp/agg.L_corp,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Capital-labor ratio in corp. sector")
if do_save==1; print(fullfile(FigDir,'Kc_Lc'),'-dpng'); end

figure
plot(agg_tran.Y_small,"linewidth",2)
yline(agg.Y_small,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Ysmall (output+liq-entry for small firms)")
if do_save==1; print(fullfile(FigDir,'Y_small'),'-dpng'); end

figure
plot(agg_tran.Y_small,"linewidth",2)
yline(agg.Y_small,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Ysmall (output+liq-entry for small firms)")
if do_save==1; print(fullfile(FigDir,'Y_small'),'-dpng'); end

%% results for eta = 1
agg_tran = results(2).agg_tran;
path     = results(2).path;

figure
plot(agg_tran.b_grid(1,:,1),'-o',"linewidth",2)
hold on
plot(agg_tran.b_grid(1,:,2),'-o',"linewidth",2)
yline(b_grid(1),'--',{'Steady-state'})
xlabel("Time in transition, t")
title("b tilde, first point of dynamic b grid")
legend("Impacted","Unimpacted")
hold off
if do_save==1; print(fullfile(FigDir,'b_grid'),'-dpng'); end

figure
plot(agg_tran.K_agg,"linewidth",2)
yline(agg.K_agg,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Aggregate capital K")
if do_save==1; print(fullfile(FigDir,'K'),'-dpng'); end

figure
plot(path.C,"linewidth",2)
yline(agg.C_agg,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Consumption")
if do_save==1; print(fullfile(FigDir,'C'),'-dpng'); end

figure
plot(path.KL_ratio,"linewidth",2)
yline(agg.K_corp/agg.L_corp,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Capital-labor ratio in corp. sector")
if do_save==1; print(fullfile(FigDir,'Kc_Lc'),'-dpng'); end

figure
plot(agg_tran.Y_small,"linewidth",2)
yline(agg.Y_small,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Ysmall (output+liq-entry for small firms)")
if do_save==1; print(fullfile(FigDir,'Y_small'),'-dpng'); end

figure
plot(agg_tran.Y_small,"linewidth",2)
yline(agg.Y_small,'--',{'Steady-state'})
xlabel("Time in transition, t")
title("Ysmall (output+liq-entry for small firms)")
if do_save==1; print(fullfile(FigDir,'Y_small'),'-dpng'); end