%% This script generates plots

% Make directory where to save the plots 
if ~(isfolder(FigDir))
    mkdir(FigDir)
end

% figure
% plot(x_grid,x_prob,'-o','linewidth',2)
% xlabel("Productivity, x")
% title("Invariant distribution of g(x'|x)")
% 
% figure
% plot(x_grid,x0_prob,'-o','linewidth',2)
% xlabel("Productivity, x")
% title("Productivity distribution of new entrants")

figure
plot(x_grid,x_prob,'-o','linewidth',2)
hold on
plot(x_grid,x0_prob,'-o','linewidth',2)
xlabel("Productivity, x")
legend('Invariant distrib','Distrib new entrants')
title("Productivity distribution")
if do_save==1; print(fullfile(FigDir,'x_prob'),'-dpng'); end

%% VFI plots
[b_grid_mat,x_grid_mat] = meshgrid(b_grid,x_grid);
% Static profit \pi(x)
figure
plot(x_grid,sol.profit_vec,'linewidth',2);
yline(0,'--','Zero');
xlim([x_grid(1),x_grid(end)])
xlabel("Productivity, x",'fontsize',14)
ylabel("Profit")
title("Static profit \pi(x)",'fontsize',14)
if do_save==1; print(fullfile(FigDir,'profit'),'-dpng'); end

figure
plot(b_grid,b_grid,'k--','linewidth',2);
hold on
plot(b_grid,sol.pol_debt(:,1),'r','linewidth',2);
hold on
plot(b_grid,sol.pol_debt(:,round(nx/2)),'g','linewidth',2);
hold on
plot(b_grid,sol.pol_debt(:,nx),'b','linewidth',2);
xlabel("Current Debt, b",'fontsize',14)
ylabel("Next-period Debt, b' ",'fontsize',14)
legend("45 line","low-x","med-x","high-x","location","best")
axis tight
title("Debt policy function b'(x,b)",'fontsize',14)
hold off
if do_save==1; print(fullfile(FigDir,'pol_debt'),'-dpng'); end
    
% % v and v0 of unconstrained firms
% figure
% plot(x_grid,sol.val_unc(10,:),'linewidth',2);
% hold on
% plot(x_grid,sol.val0_unc(10,:),'linewidth',2);
% legend('v','v0')
% yline(theta-b_grid(10),'--','Resale value \theta - b');
% xlim([x_grid(1),x_grid(end)])
% xlabel("Productivity, x",'fontsize',14)
% title("Value functions for unconstraint firms",'fontsize',14)
% hold off

% % v and v0 of all firms
% figure
% plot(x_grid,sol.val(10,:),'linewidth',2);
% hold on
% plot(x_grid,sol.val0(10,:),'linewidth',2);
% legend('v','v0')
% yline(theta-b_grid(10),'--','Resale value \theta - b');
% xlim([x_grid(1),x_grid(end)])
% xlabel("Productivity, x",'fontsize',14)
% title("Value functions for all firms",'fontsize',14)
% hold off

 % v and v0 unconstraint firms, 3-D plot
figure
surf(b_grid_mat,x_grid_mat,sol.val_unc')
hold on 
surf(b_grid_mat,x_grid_mat,sol.val0_unc')
xlabel('Current Debt, b','Fontsize',14)
ylabel('Current Productivity, x','Fontsize',14)
axis tight
zlabel('Value, v(b,x)','Fontsize',14)
title('Value function for unconstraint firms')
hold off
if do_save==1; print(fullfile(FigDir,'val_unc'),'-dpng'); end

% v and v0 all firms, 3-D plot
figure
surf(b_grid_mat,x_grid_mat,sol.val')
hold on 
surf(b_grid_mat,x_grid_mat,sol.val0')
xlabel('Current Debt, b','Fontsize',14)
ylabel('Current Productivity, x','Fontsize',14)
axis tight
zlabel('Value, v(b,x)','Fontsize',14)
title('Value function for firms')
hold off
if do_save==1; print(fullfile(FigDir,'val_all'),'-dpng'); end

% Exit policy   COLOR MAP
figure
%contourf(a_grid,theta_grid,occ,[1 2 3 4])
imagesc(b_grid,x_grid,sol.pol_exit','CDataMapping','scaled')
set(gca,'YDir','normal')  %important line
colorbar
xlabel('Debt, b','fontsize',14)
ylabel('Productivity, x','fontsize',14)
title('1 = EXIT, 0 = STAY','fontsize',14)
if do_save==1; print(fullfile(FigDir,'exit'),'-dpng'); end

% Entry policy  COLOR MAP
figure
%contourf(a_grid,theta_grid,occ,[1 2 3 4])
imagesc(b_grid,x_grid,sol.pol_entry','CDataMapping','scaled')
set(gca,'YDir','normal')  %important line
colorbar
xlabel('Debt, b','fontsize',14)
ylabel('Productivity, x','fontsize',14)
title('1 = Entry, 0 = No Entry','fontsize',14)
if do_save==1; print(fullfile(FigDir,'entry'),'-dpng'); end

%% plots in mu file

% Distribution mu_0, 3-D plot
figure
surf(b_grid_mat,x_grid_mat,mu')
xlabel('Current Debt, b','Fontsize',14)
ylabel('Current Productivity, x','Fontsize',14)
axis tight
zlabel('Mu, mu(b,x)','Fontsize',14)
title('Distribution of firms')
if do_save==1; print(fullfile(FigDir,'mu_3d'),'-dpng'); end

% Distribution mu, 2-D plot
figure
plot(b_grid,mu(:,1),'linewidth',2)
hold on
plot(b_grid,mu(:,nx),'linewidth',2)
legend('Lowest x','Highest x')
xlabel('Current Debt, b','Fontsize',14)
ylabel('mu(b,x)','Fontsize',14)
%axis tight
title('Distribution')
hold off

% Marginal distrib of debt mu_b
mu_debt = sum(mu,2);
plot(b_grid,mu_debt,'linewidth',2)
xlabel('Current Debt, b','Fontsize',14)
ylabel('mu(b)','Fontsize',14)
axis tight
title('Marginal PDF of debt b')
if do_save==1; print(fullfile(FigDir,'mu_debt'),'-dpng'); end

