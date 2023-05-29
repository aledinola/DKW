function [] = plot_ss_policy(b_grid,sol,distribS,par,SaveDir,FormatFig,do_save,FS)

% DESCRIPTION
% plot_ss_policy plots policy function of investment 
% inv(k,b,x):= k'(k,b,x)-(1-delta)k
%
% INPUTS
% - b_grid     : steady state b_grid
% - sol        : steady state policy functions, structure
% - distribS   : steady state distributions, structure
% - par        : parameters, structure
% - SaveDir    : directory for output figures
% - FormatFig  : format of output figures
% - do_save    : 0/1
% - FS         : Font size
% OUTPUTS
% plot_ss_policy generates and saves the following figures:
% - inv_type_ss
% - 

% Unpack
pol_exit = sol.pol_exit;
pol_debt = sol.pol_debt; % (nk,nb,nx)
pol_kp = sol.pol_kp;  % (nk,nb,nx)
B_hat  = sol.B_hat;   % (nk,nx)
k_grid = par.k_grid;  % (nk,1)
x_grid = par.x_grid;
delta  = par.delta_k; % scalar
mu     = distribS.mu;
[nk,nb,nx] = size(pol_kp);

inv = zeros(nk,nb,nx);
inv_ave = zeros(nk,nx); % ave investment for each (b,x)
inv_tot = zeros(nk,nx); % total investment for each (b,x)
inv_rate = zeros(nk,nx); % investment rate for each (b,x)
mu_tot = zeros(nk,nx);
inv_type = zeros(nk,nb,nx); % 0 = exit, 1 = disinvest, 2 = invest (or inaction), 3 = spike invest
borrow_type = zeros(nk,nb,nx); % 0 = firm saves, 1 = firm borrows
borrow_pol  = zeros(nk,nb,nx); % pol_debt - b
for x_c = 1:nx
    for k_c = 1:nk
        k_val = k_grid(k_c);

        for b_c = 1:nb

            b_val = b_grid(k_c,b_c);
            borrowing = pol_debt(k_c,b_c,x_c)-b_val;
            borrow_pol(k_c,b_c,x_c) = borrowing;
            %inv(k_c,b_c,x_c) = pol_kp(k_c,b_c,x_c)-(1-delta)*k_val;
            inv(k_c,b_c,x_c) = (1-pol_exit(k_c,b_c,x_c))*pol_kp(k_c,b_c,x_c)-(1-delta)*k_val;

            if pol_exit(k_c,b_c,x_c)<0.5 &&  pol_kp(k_c,b_c,x_c)-(1-delta)*k_val<0
                inv_type(k_c,b_c,x_c) = 1;
            elseif pol_exit(k_c,b_c,x_c)<0.5 &&  pol_kp(k_c,b_c,x_c)-(1-delta)*k_val>=0 ...
                    &&  (pol_kp(k_c,b_c,x_c)-(1-delta)*k_val)/k_val<1
                inv_type(k_c,b_c,x_c) = 2;
            elseif pol_exit(k_c,b_c,x_c)<0.5 &&  pol_kp(k_c,b_c,x_c)-(1-delta)*k_val>=0 ...
                    &&  (pol_kp(k_c,b_c,x_c)-(1-delta)*k_val)/k_val>=1
                inv_type(k_c,b_c,x_c) = 3;
            end
            %invrate(k_c,b_c,x_c) = pol_kp(k_c,b_c,x_c)/k_val-(1-delta);
            
            if borrowing<=0 
                borrow_type(k_c,b_c,x_c) = 0; % firm saves
            else
                borrow_type(k_c,b_c,x_c) = 1; % firm borrows
            end
        end
        inv_tot(k_c,x_c) = sum(inv(k_c,:,x_c) .* mu(k_c,:,x_c),'all');
        mu_tot(k_c,x_c) =  sum(mu(k_c,:,x_c),'all');
        inv_ave(k_c,x_c) = inv_tot(k_c,x_c) / sum(mu(k_c,:,x_c),'all');
        inv_rate(k_c,x_c) = sum(inv(k_c,:,x_c) .* mu(k_c,:,x_c),'all')/ ...
            sum(k_val .* mu(k_c,:,x_c),'all');
    end
end
%mu_k  = sum(mu_tot,2);

%% Make plot
indb = ones(nk,1);
for k_c = 1:nk
    [~,indb(k_c)] = min(abs(b_grid(k_c,:)+50));
end


% figure('Position', [1, 1, 800, 300])
% tiledlayout(1,length(20:30:50))
% for k_c = 20:30:50
%     nexttile
%     contourf(b_grid(k_c,indb(k_c)-1:end),x_grid,squeeze(inv(k_c,indb(k_c)-1:end,:))','LineStyle','none')
%     %yline(x_tilde_val(k_c),'color','w','LineWidth',2)
%     xlabel('Debt, b','fontsize',FS)
%     ylabel('Productivity, x','fontsize',FS)
%     title(['k = ',num2str(par.k_grid(k_c),'%4.2f')],'fontsize',FS)
%     %txt1 = 'Entry';
%     %txt2 = 'No Entry';
%     %text(-30,2.2,txt1,'FontSize',20,'Color', 'k','FontWeight','bold')
%     %text(-30,0.5,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
% end
% if do_save==1; print(fullfile(SaveDir,'inv_ss'),FormatFig); end

figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
    nexttile
    contourf(b_grid(k_c,indb(k_c)-1:end),x_grid,squeeze(inv_type(k_c,indb(k_c)-1:end,:))',[0 1 2  3])
    %yline(x_tilde_val(k_c),'color','w','LineWidth',2)
    xlabel('Debt, b','fontsize',FS)
    ylabel('Productivity, x','fontsize',FS)
    title(['k = ',num2str(par.k_grid(k_c),'%4.2f')],'fontsize',FS)
    %txt1 = 'Entry';
    %txt2 = 'No Entry';
    %text(-30,2.2,txt1,'FontSize',20,'Color', 'k','FontWeight','bold')
    %text(-30,0.5,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
end
if do_save==1; print(fullfile(SaveDir,'inv_type_ss'),FormatFig); end

%% Borrowing decision

% figure('Position', [1, 1, 800, 300])
% tiledlayout(1,length(20:30:50))
% for k_c = 20:30:50
%     nexttile
%     contourf(b_grid(k_c,indb(k_c)-1:end),x_grid,squeeze(borrow_type(k_c,indb(k_c)-1:end,:))',[0 1])
%     %yline(x_tilde_val(k_c),'color','w','LineWidth',2)
%     xlabel('Debt, b','fontsize',FS)
%     ylabel('Productivity, x','fontsize',FS)
%     title(['k = ',num2str(par.k_grid(k_c),'%4.2f')],'fontsize',FS)
%     %txt1 = 'Entry';
%     %txt2 = 'No Entry';
%     %text(-30,2.2,txt1,'FontSize',20,'Color', 'k','FontWeight','bold')
%     %text(-30,0.5,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
% end
% if do_save==1; print(fullfile(SaveDir,'borrow_type_ss'),FormatFig); end


figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
    nexttile
    contourf(b_grid(k_c,indb(k_c)-1:end),x_grid,squeeze(borrow_pol(k_c,indb(k_c)-1:end,:))')
    %yline(x_tilde_val(k_c),'color','w','LineWidth',2)
    xlabel('Debt, b','fontsize',FS)
    ylabel('Productivity, x','fontsize',FS)
    title(['k = ',num2str(par.k_grid(k_c),'%4.2f')],'fontsize',FS)
    %txt1 = 'Entry';
    %txt2 = 'No Entry';
    %text(-30,2.2,txt1,'FontSize',20,'Color', 'k','FontWeight','bold')
    %text(-30,0.5,txt2,'FontSize',20,'Color', 'w','FontWeight','bold')
end
if do_save==1; print(fullfile(SaveDir,'borrow_pol_ss'),FormatFig); end

%% Fraction of unconstrained firms

is_unc = zeros(nk,nb,nx);
for x_c = 1:nx
    for b_c = 1:nb
        for k_c = 1:nk
            if b_grid(k_c,b_c)<=B_hat(k_c,x_c)
                is_unc(k_c,b_c,x_c) = 1;
            end
        end
    end
end
frac_unc = sum(is_unc.*mu,'all')/sum(mu,'all');
if do_save==1
    FID = fopen(fullfile('tables','frac_unc.txt'),'w');
    fprintf(FID,'Fraction of unconstrained firms: %8.4f \n',frac_unc);
    fclose(FID);
end

end %end function "plot_ss_policy"