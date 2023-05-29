
if par.grant_flag==1 && par.grant_target == 0
    par.FigDir   = fullfile('figures','grant_baseline');
elseif par.grant_flag==5 && par.grant_target == 1
    par.FigDir   = fullfile('figures','grant_targeted');
elseif par.grant_flag==0
    par.FigDir   = fullfile('figures','nogrant');
elseif par.grant_flag==3 && par.grant_target == 0
    par.FigDir   = fullfile('figures','grant_large');
elseif par.grant_flag==4 && par.grant_target == 0
    par.FigDir   = fullfile('figures','grant_small');
end




%---------------- Shocks ----------------%
    figure
    plot(1:par.T+1,par.A_small(:,1),"-o","linewidth",2)
    hold on
    plot(1:par.T+1,par.A_corp,"-o","linewidth",2)
    hold on
    plot(1:par.T+1,par.margutil,"-o","linewidth",2)
    hold on
    plot(1:par.T+1,par.lsupply,"-o","linewidth",2)
    yline(1,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    legend("A_{small}","A_{corp}","Util","L-Sup")
    title("Shocks")
    hold off
    if do_save==1; print(fullfile(par.FigDir,'shocks'),'-dpng'); end
    


%---------------- Representative agent part of the model ----------------%
    figure
    plot(results.path.C,"linewidth",2)
    hold on
    yline(agg.C_agg,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Consumption")
    if do_save==1; print(fullfile(par.FigDir,'C'),'-dpng'); end
    
    figure
    plot(results.path.KL_ratio,"linewidth",2)
    hold on
    yline(agg.K_corp/agg.L_corp,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Capital-labor ratio in corp. sector")
    if do_save==1; print(fullfile(par.FigDir,'Kc_Lc'),'-dpng'); end
    
    figure
    plot(results.path.q,"linewidth",2)
    hold on
    yline(prices.q,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Financial discount factor, q")
    if do_save==1; print(fullfile(par.FigDir,'q'),'-dpng'); end
    
    figure
    plot(results.path.rental,"linewidth",2)
    hold on
    yline(prices.rental,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Rental rate, R")
    if do_save==1; print(fullfile(par.FigDir,'rental'),'-dpng'); end
    
    
    figure
    plot(results.path.w,"linewidth",2)
    hold on
    yline(prices.wage,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Wage")
    if do_save==1; print(fullfile(par.FigDir,'wage'),'-dpng'); end
    
    %------------------- Heterog agents --------------------------%
    % dim is (b,kappa,time,impact x grant)
    for k_c = 1:par.nk
        mytitle = ['b tilde, k=',num2str(k_c),', impacted with grant'];
        mysave = ['b_grid_k',num2str(k_c)];
        figure
        plot(squeeze(results.agg_tran.b_grid(1,k_c,:,1)),"linewidth",2)
        hold on
        yline(b_grid(1,k_c),'--',{'Steady-state'});
        xlabel("Time in transition, t")
        title(mytitle)
        hold off
        if do_save==1; print(fullfile(par.FigDir,mysave),'-dpng'); end
    end
    
    
    %%% Household Capital
    figure
    plot(results.agg_tran.K_agg,"linewidth",2)
    hold on
    yline(agg.K_agg,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Household capital")
    if do_save==1; print(fullfile(par.FigDir,'K_agg'),'-dpng'); end
    
    %%% Economy-wide Capital
    figure
    plot(results.agg_tran.K_all,"linewidth",2)
    hold on
    yline(agg.K_all,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Aggregate capital")
    if do_save==1; print(fullfile(par.FigDir,'K_all'),'-dpng'); end
    
    figure
    plot(results.agg_tran.K_corp,"linewidth",2)
    hold on
    yline(agg.K_corp,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Capital corporate sector")
    if do_save==1; print(fullfile(par.FigDir,'K_corp'),'-dpng'); end
    
    figure
    plot(results.agg_tran.K_small,"linewidth",2)
    hold on
    yline(agg.K_small,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Capital small firms")
    if do_save==1; print(fullfile(par.FigDir,'K_small'),'-dpng'); end
    
    figure
    plot(results.agg_tran.mass_small,"linewidth",2)
    hold on
    yline(agg.mass_small,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Mass of small firms")
    if do_save==1; print(fullfile(par.FigDir,'mass_small'),'-dpng'); end
    
    %%% Output
    figure
    plot(results.agg_tran.Y_agg,"linewidth",2)
    hold on
    yline(agg.Y_agg,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Aggregate Output")
    if do_save==1; print(fullfile(par.FigDir,'Y_agg'),'-dpng'); end
    
    figure
    plot(results.agg_tran.output_small,"linewidth",2)
    hold on
    yline(agg.output_small,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Output small firms")
    if do_save==1; print(fullfile(par.FigDir,'output_small'),'-dpng'); end
    
    figure
    plot(results.agg_tran.Y_corp,"linewidth",2)
    hold on
    yline(agg.Y_corp,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Output corporate sector")
    if do_save==1; print(fullfile(par.FigDir,'Y_corp'),'-dpng'); end
    
    figure
    plot(results.agg_tran.Y_small,"linewidth",2)
    hold on
    yline(agg.Y_small,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Ysmall (output+liq-entry for small firms)")
    if do_save==1; print(fullfile(par.FigDir,'Y_small'),'-dpng'); end
    
    
    
    figure
    plot(results.agg_tran.InvK_all,"linewidth",2)
    hold on
    yline(agg.InvK_all,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Capital Investment")
    if do_save==1; print(fullfile(par.FigDir,'InvK'),'-dpng'); end
    
    %%% Employment
    figure %OK
    plot(results.agg_tran.L_agg,"linewidth",2)
    hold on
    yline(agg.L_agg,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Aggregate Employment")
    if do_save==1; print(fullfile(par.FigDir,'L_agg'),'-dpng'); end
    
    figure 
    plot(results.agg_tran.L_small,"linewidth",2)
    hold on
    yline(agg.L_small,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Employment small firms")
    if do_save==1; print(fullfile(par.FigDir,'L_small'),'-dpng'); end
    
    figure 
    plot(results.agg_tran.L_corp,"linewidth",2)
    hold on
    yline(agg.L_corp,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Employment corporate sector")
    if do_save==1; print(fullfile(par.FigDir,'L_corp'),'-dpng'); end
    
    figure 
    plot(results.agg_tran.entry_vec,"linewidth",2)
    hold on
    yline(agg.entry,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Entry costs")
    if do_save==1; print(fullfile(par.FigDir,'entry'),'-dpng'); end
    
    figure
    plot(results.agg_tran.liq_vec,"linewidth",2)
    hold on
    yline(agg.liq,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Liquidation costs")
    if do_save==1; print(fullfile(par.FigDir,'liq'),'-dpng'); end

    figure
    plot(results.agg_tran.exit_rate_vec,"linewidth",2)
    hold on
    yline(agg.exit_rate,'--',{'Steady-state'});
    xlabel("Time in transition, t")
    title("Exit rate")
    if do_save==1; print(fullfile(par.FigDir,'exit'),'-dpng'); end