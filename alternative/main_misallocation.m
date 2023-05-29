%% This script generates results regarding capital misallocation
% Reference: "Misallocation.pdf"
% The script loads the following mat files:
% ss.mat, grant_baseline.mat, nogrant.mat, grant_targslim.mat
% which must be saved in subfolder "mat"
% -

clear;clc;close all
% Add path to numerical tools
addpath(genpath(fullfile('tools')));

%% Load mat files
disp('Loading distributions')
% Steady state
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'),'distribS','agg','par')
else
    error('File ss.mat does not exist!')
end
mu_ss = distribS.mu;

% Baseline grant
if isfile(fullfile('mat','grant_baseline.mat'))
    load(fullfile('mat','grant_baseline.mat'),'distrib_tran','agg_tran','par')
else
    error('File grant_baseline.mat does not exist!')
end
mu_baseline = distrib_tran.mu;
weights_baseline = par.weights;
agg_baseline = agg_tran;

% No grant
if isfile(fullfile('mat','nogrant.mat'))
    load(fullfile('mat','nogrant.mat'),'distrib_tran','agg_tran','par')
else
    error('File nogrant.mat does not exist!')
end
mu_nogrant = distrib_tran.mu;
weights_nogrant = par.weights;
agg_nogrant = agg_tran;

% Targeted grant
if isfile(fullfile('mat','grant_targslim.mat'))
    load(fullfile('mat','grant_targslim.mat'),'distrib_tran','agg_tran','par')
else
    error('File grant_targslim.mat does not exist!')
end
mu_targslim = distrib_tran.mu;
weights_targslim = par.weights;
agg_targslim = agg_tran;

%% Q1: How much TFP is lost because capital is misallocated between small firms?
disp('Q1: misallocation within small-firm sector')

% -Steady state
%[Y_small_opt_ss,fc_ss] = fun_misallocation_Q1(agg.K_small,agg.L_small,distribS.mu,par.x_grid,par);
% Compute tfp_small (optimal effective small-firm sector TFP)
nu_x = squeeze(sum(mu_ss,[1,2]));
tfp_small_ss = sum(par.x_grid.^(1/(1-par.gamma2)).*nu_x)^(1-par.gamma2);
% Compute optimal small firm output
Y_small_opt_ss = tfp_small_ss*par.A*agg.K_small^(par.gamma1*par.gamma2)*agg.L_small^(par.gamma2*(1-par.gamma1));
% Compute the actual small firm output
% Note that agg.output_small is the total output minus total fixed cost.
% We add back the total fixed cost.
fc_ss = 0;
for k_c = 1:par.nk
    fc_ss = fc_ss + par.fixcost(k_c)*sum(mu_ss(k_c,:,:),'all');
end
Y_small_actual_ss = agg.output_small + fc_ss;
% Misallocation rate in the ss
mis_rate_ss = (Y_small_opt_ss - Y_small_actual_ss)/Y_small_opt_ss;

% Baseline grant
[mis_rate_baseline,Y_small_opt_baseline,Y_small_actual_baseline,tfp_small_baseline] = ...
    fun_misallocation_Q1(agg_baseline,mu_baseline,weights_baseline,par);

% No grant
[mis_rate_nogrant,Y_small_opt_nogrant,Y_small_actual_nogrant,tfp_small_nogrant] = ...
    fun_misallocation_Q1(agg_nogrant,mu_nogrant,weights_nogrant,par);

% Targeted grant
[mis_rate_targslim,Y_small_opt_targslim,Y_small_actual_targslim,tfp_small_targslim] = ...
    fun_misallocation_Q1(agg_targslim,mu_targslim,weights_targslim,par);

%% Plot
SaveDir = fullfile('figures','grant_vs_nogrant'); %folder where figures are stored
% Options for plots
FormatFig= '-dpng'; % Specify '-dpng' or '-depsc'
do_save  = 1;        % Flag 0-1
lastp    = 40;       % number of periods plotted
FS       = 18;
LW       = 4;

% Plot (baseline vs no grant)
figure
plot(0:lastp-1,mis_rate_baseline(1:lastp),'-',"linewidth",LW)
hold on
plot(0:lastp-1,mis_rate_nogrant(1:lastp),'--o',"linewidth",LW-1)
yline(mis_rate_ss,'--',"linewidth",LW-1);
legend('Baseline grant','Laissez-faire','Steady state','FontSize',FS,'Location','Best')
ax = gca;
ax.FontSize = FS;
xlabel("Time in transition, t",'FontSize',FS)
ylabel("Misallocation rate",'FontSize',FS)
grid on
%title(label,'FontSize',myfontsize)
if do_save==1; print(fullfile(SaveDir,'mis_rate'),FormatFig); end

% Plot (baseline vs targeted)
figure
plot(0:lastp-1,mis_rate_baseline(1:lastp),'-',"linewidth",LW)
hold on
plot(0:lastp-1,mis_rate_targslim(1:lastp),'--x',"linewidth",LW-1)
yline(mis_rate_ss,'--',"linewidth",LW-1);
legend('Baseline grant','Targeted grant','Steady state','FontSize',FS,'Location','Best')
ax = gca;
ax.FontSize = FS;
xlabel("Time in transition, t",'FontSize',FS)
ylabel("Misallocation rate",'FontSize',FS)
grid on
%title(label,'FontSize',myfontsize)
if do_save==1; print(fullfile(SaveDir,'mis_rate_2'),FormatFig); end

%% Q2: How much TFP is lost because capital is misallocated between the two sectors?
disp('Q2: misallocation between sectors')

% Solve maximization problem in Q2 numerically
alpha   = par.alpha;
gamma1  = par.gamma1;
gamma2  = par.gamma2;
fixcost2 = 0;%par.fixcost2;
K_agg  = agg.K_agg;
L_agg  = agg.L_agg;
A_c    = 1; %in the steady-state
xbar   = tfp_small_ss;%*(1-0.11);

myfun = @(x) -fun_opt_Y(x,K_agg,L_agg,A_c,xbar,alpha,gamma1,gamma2,fixcost2);

x_init = [agg.K_corp,agg.L_corp]';
LB = [0,0]';
UB = [K_agg,L_agg]'; 

% ng = 100;
% x1_vec = linspace(0,K_agg,ng)';
% x2_vec = linspace(0,L_agg,ng)';
% y_mat = zeros(ng,ng);
% for i2=1:ng
%     for i1=1:ng
%         xvec = [x1_vec(i1),x2_vec(i2)]';
%         y_mat(i1,i2) = fun_opt_Y(xvec,K_agg,L_agg,A_c,xbar,alpha,gamma1,gamma2,fixcost2);
%     end
% end
% 
% [f_opt_grid,linind] = max(y_mat(:));
% [x_opt_grid_1,x_opt_grid_2] = ind2sub([ng,ng],linind);
% x_opt1 =[x1_vec(x_opt_grid_1),x2_vec(x_opt_grid_2)]';

[x_opt,f_opt] = fminsearchcon(myfun,x_init,LB,UB);

% Steady state: L_corp_guess is negative!
L_corp_guess = agg.L_agg-(tfp_small_ss*par.gamma2)^(1/(1-par.gamma2))*(agg.K_agg/agg.L_agg)^(-par.alpha);


