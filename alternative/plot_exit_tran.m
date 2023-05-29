function [exit_rule_ss,exit_rule_tran] = plot_exit_tran(par,b_grid_ss,sol,pol_tran,SaveDir,do_save,FormatFig)

% Purpose: plot exit thresholds in impact period of transition, t=1
% INPUTS:
%   "b_grid_ss" Grid for debt in the steady-state, dim: (nk,nb)
%   "sol"       Structure with steady-state policy and value functions
%   "pol_tran"  Struct with pol_exit and b_grid for transition
% OUTPUTS:
%   "exit_rule_ss"    Exit threshold steady-state, dim: (nk,nb)
%   "exit_rule_tran"  Exit threshold transition,   dim: (nk,nb)


pol_exit_ss = sol.pol_exit; %dim: (k,b,x) 
% Take only t=1
pol_exit    = squeeze(pol_tran.pol_exit(:,:,:,1,:)); % dim: (k,b,x,impact x grant)
%mu_active  = squeeze(distrib_tran.mu_active(:,:,1,:));
b_grid      = squeeze(pol_tran.b_grid(:,:,1,:)); % dim: (k,b,impact x grant)

%% exit_rule: threshold x (on grid) such that firm will not exit.
% steady state
exit_rule_ss = zeros(par.nk,par.nb);
for b_c = 1:par.nb
    for k_c = 1:par.nk
        % Smallest x s.t. firm does not exit
        ind = find(squeeze(pol_exit_ss(k_c,b_c,:))==0,1,'first');
        if ~isempty(ind)
            exit_rule_ss(k_c,b_c) = par.x_grid(ind);
        else
            exit_rule_ss(k_c,b_c) = par.x_grid(par.nx)+1e-10;
        end
        % If ind is empty, out of bounds error
        
    end
end
% transition
exit_rule_tran = zeros(par.nk,par.nb,par.nn);
for n_c = 1:par.nn % impact * grant
    for b_c = 1:par.nb
        for k_c = 1:par.nk
            ind = find(squeeze(pol_exit(k_c,b_c,:,n_c))==0,1,'first');
            if ~isempty(ind)
                exit_rule_tran(k_c,b_c,n_c) = par.x_grid(ind);
            else
                exit_rule_tran(k_c,b_c,n_c) = par.x_grid(par.nx)+1e-10;
            end
        end
    end
end

figure
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
% 1 = impact, grant,  2 = no impact, grant, 3 = impact, no grant, 4 =  no impact, no grant
% b_grid is (k,b,n)
% exit_rule_tran is (k,b,n)
    nexttile
    plot(squeeze(b_grid(k_c,:,1)),squeeze(exit_rule_tran(k_c,:,1)),"linewidth",3);
    hold on
    plot(squeeze(b_grid(k_c,:,3)),squeeze(exit_rule_tran(k_c,:,3)),'--',"linewidth",3);
    hold on 
    plot(b_grid_ss(k_c,:),exit_rule_ss(k_c,:),"linewidth",3);
    hold on
    legend({'Impacted with grant','Impacted no grant','Steady-state'},'FontSize',14,'Location','Best');
    xlabel("Debt, b")
    ylabel('Productivity, x')
    xlim([b_grid(k_c,1),b_grid(k_c,end)])
    ylim([par.x_grid(1),par.x_grid(end)])
    grid on
    % txt = '\downarrow Exit threshold';
    % text(-30,2.2,txt,'FontSize',14)
    hold off
end
if do_save==1; print(fullfile(SaveDir,'exit_tran_impacted'),FormatFig); end

figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
% 1 = impact, grant,  2 = no impact, grant, 3 = impact, no grant, 4 =  no impact, no grant
    nexttile
    plot(squeeze(b_grid(k_c,:,2)),squeeze(exit_rule_tran(k_c,:,2)),"linewidth",3);
    hold on
    plot(squeeze(b_grid(k_c,:,4)),squeeze(exit_rule_tran(k_c,:,4)),'--',"linewidth",3);
    hold on 
    plot(b_grid_ss(k_c,:),exit_rule_ss(k_c,:),"linewidth",3);
    hold on
    legend({'Not impacted with grant','Not impacted, no grant','Steady-state'},'FontSize',14,'Location','Best');
    xlabel("Debt, b")
    ylabel('Productivity, x')
    xlim([b_grid(k_c,1),b_grid(k_c,end)])
    ylim([par.x_grid(1),par.x_grid(end)])
    grid on
    % txt = '\downarrow Exit threshold';
    % text(-30,2.2,txt,'FontSize',14)
    title(['\kappa = ',num2str(par.k_grid(k_c))],'fontsize',14)
    hold off
end
if do_save==1; print(fullfile(SaveDir,'exit_tran_unimpacted'),FormatFig); end

% We do the same, but zoomed in at positive values for b
indb = ones(par.nk,1);
for k_c = 1:par.nk
    [~,indb(k_c)] = min(abs(b_grid_ss(k_c,:)-0));
end

figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
% 1 = impact, grant,  2 = no impact, grant, 3 = impact, no grant, 4 =  no impact, no grant
% b_grid is (k,b,n)
% exit_rule_tran is (k,b,n)
    nexttile
    plot(squeeze(b_grid(k_c,indb(k_c)-1:end,1)),squeeze(exit_rule_tran(k_c,indb(k_c)-1:end,1)),"linewidth",3);
    hold on
    plot(squeeze(b_grid(k_c,indb(k_c)-1:end,3)),squeeze(exit_rule_tran(k_c,indb(k_c)-1:end,3)),'--',"linewidth",3);
    hold on 
    plot(b_grid_ss(k_c,indb(k_c)-1:end),exit_rule_ss(k_c,indb(k_c)-1:end),"linewidth",3);
    hold on
    legend({'Impacted with grant','Impacted no grant','Steady-state'},'FontSize',14,'Location','Best');
    xlabel("Debt, b")
    ylabel('Productivity, x')
    xlim([b_grid(k_c,indb(k_c)),b_grid(k_c,end)])
    ylim([par.x_grid(1),par.x_grid(end)])
    grid on
    % txt = '\downarrow Exit threshold';
    % text(-30,2.2,txt,'FontSize',14)
    title(['\kappa = ',num2str(par.k_grid(k_c))],'fontsize',14)

    hold off
end
if do_save==1; print(fullfile(SaveDir,'exit_tran_cut_impacted'),FormatFig); end

figure('Position', [1, 1, 800, 300])
tiledlayout(1,length(20:30:50))
for k_c = 20:30:50
% 1 = impact, grant,  2 = no impact, grant, 3 = impact, no grant, 4 =  no impact, no grant
    nexttile
    plot(squeeze(b_grid(k_c,indb(k_c)-1:end,2)),squeeze(exit_rule_tran(k_c,indb(k_c)-1:end,2)),"linewidth",3);
    hold on
    plot(squeeze(b_grid(k_c,indb(k_c)-1:end,4)),squeeze(exit_rule_tran(k_c,indb(k_c)-1:end,4)),'--',"linewidth",3);
    hold on 
    plot(b_grid_ss(k_c,indb(k_c)-1:end),exit_rule_ss(k_c,indb(k_c)-1:end),"linewidth",3);
    hold on
    legend({'Not impacted with grant','Not impacted, no grant','Steady-state'},'FontSize',14,'Location','Best');
    xlabel("Debt, b")
    ylabel('Productivity, x')
    xlim([b_grid(k_c,indb(k_c)),b_grid(k_c,end)])
    ylim([par.x_grid(1),par.x_grid(end)])
    grid on
    % txt = '\downarrow Exit threshold';
    % text(-30,2.2,txt,'FontSize',14)
    title(['\kappa = ',num2str(par.k_grid(k_c))],'fontsize',14)

    hold off
end
if do_save==1; print(fullfile(SaveDir,'exit_tran_cut_unimpacted'),FormatFig); end

end % end function <plot_exit_tran>
