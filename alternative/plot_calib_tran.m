function [] = plot_calib_tran(model_mom_trans,par,FORMAT,FS,do_save,SaveDir)

T_last = 4; % plot the first four quarters
%% Plot shocks



eta_i         = par.eta_i; % fraction of impacted small firms
v_corp        = par.v_corp; % TFP shock corporate sector
util_shift    = par.util_shift; % Demand shock
lsupply_shift = par.lsupply_shift; % Labor supply shock
rho_shock     = par.rho_shock; % Autocorrelation of shocks
v_small       = par.v_small;
lambda_shift  = par.lambda_shift;
mass_shift    = par.mass_shift;

v_corp_vec = zeros(T_last,1);
v_small_vec = zeros(T_last,1);
util_shift_vec = zeros(T_last,1);
lsupply_shift_vec = zeros(T_last,1);
lambda_shift_vec = zeros(T_last,1);
mass_shift_vec = zeros(T_last,1);
for t_c = 1:T_last
    v_corp_vec(t_c) = v_corp*rho_shock^(t_c-1);
    v_small_vec(t_c) = v_small*rho_shock^(t_c-1);
    util_shift_vec(t_c) = util_shift*rho_shock^(t_c-1);
    lsupply_shift_vec(t_c) = lsupply_shift*rho_shock^(t_c-1);
    lambda_shift_vec(t_c) = lambda_shift*rho_shock^(t_c-1);
    mass_shift_vec(t_c) = mass_shift*rho_shock^(t_c-1);
end


figure('Position',[1000 800 560 550])
plot(1:T_last,v_corp_vec,'-','linewidth',3)
hold on
plot(1:T_last,v_small_vec,'--','linewidth',3)
hold on
plot(1:T_last,lambda_shift_vec,':','linewidth',3)
hold on
plot(1:T_last,mass_shift_vec,'-.','linewidth',3)
hold on
plot(1:T_last,util_shift_vec,'-','linewidth',3)
hold on
plot(1:T_last,lsupply_shift_vec,'--','linewidth',3)
hold on
xlabel('Quarter','FontSize',FS)
xticks(1:T_last)
xticklabels({'2020Q2','2020Q3','2020Q4','2021Q1'})
%title('Effect of grant on exit rate','FontSize',14)
legend('TFP shock, corp. \nu^c','TFP shock, impacted small firms \nu^n', ...
    'Credit shock to small firms \nu^\lambda','Shock to mass of entrants', ...
    'Demand shock \nu^d','Labor supply shock \nu^l','FontSize',FS-1,'location','southoutside')
hold off
if do_save==1; print(fullfile(SaveDir,'calib_shocks'),FORMAT); end
% 

%% Plot targeted transition moments


figure('Position',[1000 800 560 550])
yyaxis left
plot(1:T_last,model_mom_trans(1,:),'-o','linewidth',2,'MarkerSize',12) % Total output
hold on
plot(1:T_last,model_mom_trans(2,:),'--o','linewidth',2,'MarkerSize',12) % Consumption
hold on
plot(1:T_last,model_mom_trans(5,:),':o','linewidth',2,'MarkerSize',12) % total employment
hold on
plot(1:T_last,model_mom_trans(6,:),'-.o','linewidth',2,'MarkerSize',12) %  employment small firms
hold on
plot(1:T_last,model_mom_trans(10,:),'-x','linewidth',2,'MarkerSize',12) %  employment entry rate
hold on
ylim([-18 18])
ylabel('% Change','FontSize',FS)
yyaxis right
plot(1:T_last,model_mom_trans(8,:),'--x','linewidth',2,'MarkerSize',12) %  employment exit rate
ylim([-40 40])
yline(0,'--');
xlabel('Quarter','FontSize',FS)
ylabel('% Change','FontSize',FS)
xticks(1:T_last)
xticklabels({'2020Q2','2020Q3','2020Q4','2021Q1'})
%title('Effect of grant on exit rate','FontSize',14)
legend('Output','Consumption', 'Employment','Employment in small firms',...
    'Small firm entry rate','Small firm exit rate (right axis)', ...
    'FontSize',FS-1,'location','southoutside')
hold off
if do_save==1; print(fullfile(SaveDir,'model_mom_trans'),FORMAT); end

end