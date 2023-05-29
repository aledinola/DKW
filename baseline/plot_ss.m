function [] = plot_ss(b_grid,sol,par,distribS,model_mom,data_mom,SaveDir,do_save,FormatFig,FS)

% DESCRIPTION
% We plot steady-state distributions, entry and exit policies, employment
% bins
%
% INPUTS
% b_grid    :  Flexible grid for debt in s.s., (nk,nb)
% sol       :  Struct with s.s. policy and value functions
% par       :  Struct with parameters
% distribS  : Struct with distributions, fields: mu, mu_active, (nk,nb,nx)
% model_mom : Struct with model moments
% data_mom  : Struct with data moments
% SaveDir   : Directory where to save plots, character
% do_save   : 0-1 flag, scalar
% FormatFig : Format of plots (png or eps), character
% FS        : Font size

% Unpack objects
nk        = par.nk; % No. of grid points for fixed capital small firms
nx        = par.nx;
nb        = par.nb;
k_grid    = par.k_grid;
mu_active = distribS.mu_active;
pol_debt  = sol.pol_debt;

TabDir = 'tables';

%pt = 0;      % Flag: 0=don't show title, 1=show title

%% Debt to asset ratio distribution
% 1. Compute debt to asset ratio for each k,b,x. (dkratio with dimension nk,nb,nx)
% 2. convert dkratio and mu into vectors, plot bars.
% 3. cap dkratio between -2 and 3
% What fraction of firms has savings greater than 3 times asset?

dk_ratio = zeros(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            %dk_ratio(k_c,b_c,x_c) = b_grid(k_c,b_c)/k_grid(k_c);
            dk_ratio(k_c,b_c,x_c) = pol_debt(k_c,b_c,x_c)/k_grid(k_c);
        end
    end
end

% Consider only firms with zero or positive debt
dk_ratio_pos = dk_ratio(4:end,:,1:end);
mu_active_pos = mu_active(4:end,:,1:end);
dk_ratio_pos = dk_ratio_pos(pol_debt(4:end,:,1:end)>=0);
mu_active_pos = mu_active_pos(pol_debt(4:end,:,1:end)>=0);

dk_ratio_largek = dk_ratio(17:end,:,:);
mu_active_largek = mu_active(17:end,:,:);

dk_ratio_neg = dk_ratio(4:end,:,1:end);
mu_active_neg = mu_active(4:end,:,1:end);
dk_ratio_neg = dk_ratio_neg(pol_debt(4:end,:,1:end)<0);
mu_active_neg = mu_active_neg(pol_debt(4:end,:,1:end)<0);

% vectorize dkratio and mu_active, sort both vectors according to dkratio.
% Compute the cumulative distribution of dkratio.
% then we can compute some percentiles of the dkratio distribution

% Required percentiles of the distribution
perc_vec       = [0.10,0.25,0.50,0.75,0.90,0.95]';
% Unconditional distribution of debt-to-assets ratio
dk_ratio_stats = quantili(dk_ratio(:),mu_active(:),perc_vec);
dk = [perc_vec,dk_ratio_stats];
disp(dk)
writematrix(dk,fullfile(TabDir,'dk_ratio_distrib.txt'),'Delimiter','tab')
% Conditional distribution given debt is positive
dk_ratio_pos_stats = quantili(dk_ratio_pos(:),mu_active_pos(:),perc_vec);
dk_pos = [perc_vec,dk_ratio_pos_stats];
disp(dk_pos)
writematrix(dk_pos,fullfile(TabDir,'dk_ratio_pos_distrib.txt'),'Delimiter','tab')
% Conditional distribution given capital is large 
dk_ratio_largek_stats = quantili(dk_ratio_largek(:),mu_active_largek(:),perc_vec);
dk_largek = [perc_vec,dk_ratio_largek_stats];
disp(dk_largek)
writematrix(dk_largek,fullfile(TabDir,'dk_ratio_largek_distrib.txt'),'Delimiter','tab')
% Conditional distribution given neg debt and large k 
dk_ratio_neg_stats = quantili(dk_ratio_neg(:),mu_active_neg(:),perc_vec);
dk_neg = [perc_vec,dk_ratio_neg_stats];
disp(dk_neg)
writematrix(dk_neg,fullfile(TabDir,'dk_ratio_neg_distrib.txt'),'Delimiter','tab')
%-------------------------------------------------------------------------%
% Distribution of entrants
% pol_entry, phi_dist have dim: (nk,nb,nx)

%% Plot exit and entry policies

% Exit policy   
% Note: bgrid    has dim: (nk,nb)
%       pol_exit has dim: (nk,nb,nx)

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
    title(['k = ',num2str(par.k_grid(k_c),'%4.2f')],'fontsize',FS)
end
% txt1 = 'Stay';
% txt2 = 'Exit';
% text(3,4.5,txt1,'FontSize',20,'Color', 'w','FontWeight','bold')
% text(3,0.4,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
if do_save==1; print(fullfile(SaveDir,'exit_type_ss'),FormatFig); end

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
    title(['k = ',num2str(par.k_grid(k_c),'%4.2f')],'fontsize',FS)
end
% txt1 = 'Stay';
% txt2 = 'Exit';
% text(3,4.5,txt1,'FontSize',20,'Color', 'w','FontWeight','bold')
% text(3,0.4,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
if do_save==1; print(fullfile(SaveDir,'entry_type_ss'),FormatFig); end

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