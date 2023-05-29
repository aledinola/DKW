lastp = min(60,par.T);


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

%---------------- Representative agent part of the model ----------------%
    figure
    plot(irf.C_agg(1:lastp)*100,"linewidth",2)
    hold on
    yline(0,'--');
    xlabel("Time in transition, t")
    ylabel('Consumption, %% change')
    if do_save==1; print(fullfile(par.FigDir,'C_irf'),'-dpng'); end
    
    figure
    plot(irf.Y_agg(1:lastp)*100,"linewidth",2)
    hold on
    yline(0,'--');
    xlabel("Time in transition, t")
    ylabel('Total output, %% change')
    if do_save==1; print(fullfile(par.FigDir,'Y_agg_irf'),'-dpng'); end
    
    figure
    plot(irf.InvK(1:lastp)*100,"linewidth",2)
    hold on
    yline(0,'--');
    xlabel("Time in transition, t")
    ylabel('Investment, %% change')
    if do_save==1; print(fullfile(par.FigDir,'InvK_irf'),'-dpng'); end
    
    figure
    plot(irf.output_small(1:lastp)*100,"linewidth",2)
    hold on
    yline(0,'--');
    xlabel("Time in transition, t")
    ylabel('Small Firms Output, %% change')
    if do_save==1; print(fullfile(par.FigDir,'output_small_irf'),'-dpng'); end
   