function [] = plot_ss(b_grid,sol,par,prices,distribS,model_mom,data_mom,SaveDir,do_save,FormatFig,FS)

% DESCRIPTION
% We plot steady-state distributions, entry and exit policies, employment
% bins
%
% INPUTS
% b_grid    :  Flexible grid for debt in s.s., (nk,nb)
% sol       :  Struct with s.s. policy and value functions
% par       :  Struct with parameters
% prices    :  Struct with model prices, fields: w,q
% distribS  : Struct with distributions, fields: mu, mu_active, (nk,nb,nx)
% model_mom : Struct with model moments
% data_mom  : Struct with data moments
% SaveDir   : Directory where to save plots, character
% do_save   : 0-1 flag, scalar
% FormatFig : Format of plots (png or eps), character
% FS        : Font size

nk = par.nk; % No. of grid points for fixed capital small firms
nx = par.nx;
nb = par.nb;
%pt = 0;      % Flag: 0=don't show title, 1=show title

% Marginal distribution of debt by capital (integrate over x)
% mu_active, dim: (nk,nb,nx)
mu_b = sum(distribS.mu_active,3); % marginal distrib debt, dim:(nk,nb)
mu_b_norm = mu_b./sum(mu_b,2); % normalized

% Marginal distribution of capital
mu_k = sum(distribS.mu_active,[2,3]); % marginal distrib capital, dim:(nk,1)
mu_k = mu_k./sum(mu_k); % normalized

% Marginal distribution of x after exit; mu_active is mu in the draft
mu_x = squeeze(sum(distribS.mu_active,[1 2])); % marginal distrib x
mu_x = mu_x/sum(mu_x);

% Marginal distribution of x before exit; mu is mu^0 in the draft
mu_x0 = squeeze(sum(distribS.mu,[1 2]));
mu_x0 = mu_x0/sum(mu_x0);

% Joint distributions
mu_kx = squeeze(sum(distribS.mu,2)); % (nk,nx)
mu_active_kx = squeeze(sum(distribS.mu_active,2)); % (nk,nx)


% Distribution of entrants
% pol_entry, phi_dist have dim: (nk,nb,nx)
x_entry = squeeze(par.mass*sum(sol.pol_entry.*sol.phi_dist,[1 2])); % dim: (nx)

%% Plot distributions


figure
plot(par.k_grid,mu_kx(:,20)/sum(mu_kx(:,20)),'-o','linewidth',2)
hold on
plot(par.k_grid,mu_kx(:,25)/sum(mu_kx(:,25)),'-d','linewidth',2)
hold on
plot(par.k_grid,mu_kx(:,30)/sum(mu_kx(:,30)),'-*','linewidth',2)
hold off
xlabel('Capital, k','fontsize',FS)
legend('low x', 'med x','high x')
if do_save==1; print(fullfile(SaveDir,'mu_kx_ss'),FormatFig); end

figure
% x_entry has dim: (nx,1)
plot(par.x_grid,x_entry,'linewidth',2)
xlim([par.x_grid(1),par.x_grid(end)])
xlabel('Productivity, x')
ylabel('Probability mass')
title('Productivity distrib. of entrants')
if do_save==1; print(fullfile(SaveDir,'x_entry_ss'),FormatFig); end

figure
plot(par.x_grid,mu_x0,'linewidth',2)
hold on
plot(par.x_grid,mu_x,'linewidth',2)
legend('Before exit','After exit')
xlim([par.x_grid(1),par.x_grid(end)])
xlabel('Productivity, x')
ylabel('Probability mass')
title('Distribution of x')
if do_save==1; print(fullfile(SaveDir,'mu_x_ss'),FormatFig); end

figure
plot(par.k_grid,mu_k,'linewidth',2)
xlim([par.k_grid(1),par.k_grid(end)])
xlabel('Capital, k')
ylabel('Probability mass')
title('Marginal distribution of capital')
if do_save==1; print(fullfile(SaveDir,'mu_k_ss'),FormatFig); end


% Note: bgrid has dim: nk,nb
%       mu_b  has dim: nk,nb 
figure
for k_c = 10:20:nk
    plot(b_grid(k_c,:),mu_b(k_c,:),'linewidth',2)
    hold on
end
legend('\kappa low','\kappa med low','\kappa med','\kappa med high','\kappa high','location','best')
%xlim([b_grid(1),b_grid(end)])
xlabel('Debt')
ylabel('Probability mass')
title('Distribution of debt')
hold off
if do_save==1; print(fullfile(SaveDir,'mu_b_ss'),FormatFig); end

figure
for k_c = 10:20:nk
    plot(b_grid(k_c,:),mu_b_norm(k_c,:),'linewidth',2)
    hold on
end
legend('\kappa low','\kappa med low','\kappa med','\kappa med high','\kappa high','location','best')
%xlim([b_grid(1),b_grid(end)])
xlabel('Debt')
ylabel('Probability mass')
title('Distribution of debt, normalized')
hold off
if do_save==1; print(fullfile(SaveDir,'mu_b_ss_norm'),FormatFig); end

% 3D plot: joint distribution of x and k

figure
surf(par.k_grid,par.x_grid,mu_kx')
xlabel('Capital, k','fontsize',14)
ylabel('Productivity, x','fontsize',FS)
title('Mass of incumbent firms','fontsize',FS)
if do_save==1; print(fullfile(SaveDir,'mu_kx_3d'),FormatFig); end

figure
surf(par.k_grid,par.x_grid,mu_active_kx')
xlabel('Capital, k','fontsize',14)
ylabel('Productivity, x','fontsize',FS)
title('Mass of active firms','fontsize',FS)
if do_save==1; print(fullfile(SaveDir,'mu_active_kx'),FormatFig); end

%% Plot exit and entry policies

% Exit policy   
% Note: bgrid    has dim: (nk,nb)
%       pol_exit has dim: (nk,nb,nx)
figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))


for k_c = 20:30:50
    nexttile
    contourf(squeeze(b_grid(k_c,:)),par.x_grid,squeeze(sol.pol_exit(k_c,:,:))','LineStyle','none')
    xlabel('Debt, b','fontsize',FS)
    ylabel('Productivity, x','fontsize',FS)
    title(['\kappa = ',num2str(par.k_grid(k_c))],'fontsize',FS)
end
if do_save==1; print(fullfile(SaveDir,'exit_ss'),FormatFig); end



% exit_type: type of firm exit. Dim nk,nb,nx
% 0 = stay, 1 = vol liq, 2 = forced liq, 3 = both liq
exit_type = zeros(par.nk,par.nb,par.nx);
exit_type(sol.pol_exit_vol>sol.pol_exit_forced) = 1;
exit_type(sol.pol_exit_forced>sol.pol_exit_vol) = 2;
exit_type(sol.pol_exit==0) = 0;
exit_type = round(exit_type);




indb = ones(nk,1);
for k_c = 1:nk
    [~,indb(k_c)] = min(abs(b_grid(k_c,:)-0));
end


%% Exit or liquidation policy
figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
    nexttile
    contourf(b_grid(k_c,indb(k_c)-1:end),par.x_grid,squeeze(exit_type(k_c,indb(k_c)-1:end,:))',[-0.5 0.5 1.5],'LineStyle','none')
    xlabel('Debt, b','fontsize',FS)
    ylabel('Productivity, x','fontsize',FS)
    title(['\kappa = ',num2str(par.k_grid(k_c))],'fontsize',FS)
end
% txt1 = 'Stay';
% txt2 = 'Exit';
% text(3,4.5,txt1,'FontSize',20,'Color', 'w','FontWeight','bold')
% text(3,0.4,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
if do_save==1; print(fullfile(SaveDir,'exit_type_ss'),FormatFig); end

figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
    nexttile
    contourf(b_grid(k_c,indb(k_c)-1:end),par.x_grid,squeeze(distribS.mu(k_c,indb(k_c)-1:end,:))','LineStyle','none')
    xlabel('Debt, b','fontsize',FS)
    ylabel('Productivity, x','fontsize',FS)
    title(['mu^0,\kappa = ',num2str(par.k_grid(k_c))],'fontsize',FS)
end
% txt1 = 'Stay';
% txt2 = 'Exit';
% text(3,4.5,txt1,'FontSize',20,'Color', 'w','FontWeight','bold')
% text(3,0.4,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
if do_save==1; print(fullfile(SaveDir,'mu_exit_type_ss'),FormatFig); end


% Entry policy  

% entry_type: dim(nk,nb,nx)
entry_type = zeros(nk,nb,nx);
entry_type(sol.pol_entry>0.5) = 1;

% Entry policy plots
figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
    nexttile
    contourf(b_grid(k_c,indb(k_c)-1:end),par.x_grid,squeeze(entry_type(k_c,indb(k_c)-1:end,:))',[-0.5 0.5],'LineStyle','none')
    xlabel('Debt, b','fontsize',FS)
    ylabel('Productivity, x','fontsize',FS)
    title(['\kappa = ',num2str(par.k_grid(k_c))],'fontsize',FS)
end
% txt1 = 'Stay';
% txt2 = 'Exit';
% text(3,4.5,txt1,'FontSize',20,'Color', 'w','FontWeight','bold')
% text(3,0.4,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
if do_save==1; print(fullfile(SaveDir,'entry_type_ss'),FormatFig); end


% The following is obsolete: entry policy if initial debt is the same for all firms
% Find the entry policy at b0 (initial debt) for each x and k
% entry_kx = zeros(nk,nx);
% for k_c = 1:nk
%     [b0_c,omega] = find_loc(b_grid(k_c,:),par.b0);
%     for x_c = 1:nx
%         entry_kx(k_c,x_c) = omega*sol.pol_entry(k_c,b0_c,x_c)+(1-omega)*sol.pol_entry(k_c,b0_c+1,x_c);
%     end
% end
% 
% % entry_type: dim(nk,nx)
% entry_type = zeros(nk,nx);
% entry_type(entry_kx>0.5) = 1;
% 
% figure
% contourf(par.k_grid,par.x_grid,entry_kx','LineStyle','none')
% xlabel('Capital, k','fontsize',14)
% ylabel('Productivity, x','fontsize',14)
% %txt1 = 'Entry';
% %txt2 = 'No Entry';
% %text(-30,2.2,txt1,'FontSize',20,'Color', 'k','FontWeight','bold')
% %text(-30,0.5,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
% if do_save==1; print(fullfile(SaveDir,'entry_ss'),FormatFig); end
% 
% figure
% contourf(par.k_grid,par.x_grid,entry_type',[-0.5 0.5],'LineStyle','none')
% xlabel('Capital, k','fontsize',14)
% ylabel('Productivity, x','fontsize',14)
% %txt1 = 'Entry';
% %txt2 = 'No Entry';
% %text(-30,2.2,txt1,'FontSize',20,'Color', 'k','FontWeight','bold')
% %text(-30,0.5,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
% if do_save==1; print(fullfile(SaveDir,'entry_type_ss'),FormatFig); end


%%  Model Fit
% Firm size distribution
Names = {'firmshare_0_9','firmshare_10_19','firmshare_20_99','firmshare_100_499'};
firm_size_mod = struct2vec(model_mom,Names);
firm_size_dat = struct2vec(data_mom,Names);

bin_names = {'0-9';'10-19';'20-99';'100+'};

figure
b = bar([firm_size_mod,firm_size_dat],1,'grouped');
b(1).FaceColor = [0 0 0]; %black
b(2).FaceColor = [0.8 0.8 0.8]; %grey
legend('Model','Data','location','northeast','FontSize',12)
xlabel('Size Bins')
%xlim([0 20])
%ylim([-10 10])
%title('Firm size distribution')
set(gca,'FontSize',FS);
set(gca,'xticklabel',bin_names)
if do_save==1; print(fullfile(SaveDir,'firmSize'),FormatFig); end

% Emp share by firm size bins - Model Fit
Names = {'empshare_0_9','empshare_10_19','empshare_20_99','empshare_100_499'};
emp_size_mod = struct2vec(model_mom,Names);
emp_size_dat = struct2vec(data_mom,Names);

figure
b = bar([emp_size_mod,emp_size_dat],1,'grouped');
b(1).FaceColor = [0 0 0]; %black
b(2).FaceColor = [0.8 0.8 0.8]; %grey
legend('Model','Data','location','northeast','FontSize',12)
xlabel('Size Bins')
%xlim([0 20])
%ylim([-10 10])
%title('Employment share by firm size')
set(gca,'FontSize',FS);
set(gca,'xticklabel',bin_names)
if do_save==1; print(fullfile(SaveDir,'empShare'),FormatFig); end

end % end function <plot_ss>