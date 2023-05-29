function [emp_decile_x,output_decile_x,ave_change_exit_x,share_change_exit_x] =...
    fun_deciles_x(nbins,exitSS,exit_nogrant,exit_grant_baseline,par)

pdf_x = exitSS.mass_x/sum(exitSS.mass_x); % pdf of x in the steady state
cum_x = cumsum(pdf_x);
%d_vec   = min(nbins,max(1,ceil(cum_x*nbins))); % dim (nx, 1)

%xcum = 0;
emp_decile_x = zeros(nbins,1);
output_decile_x = zeros(nbins,1);
ave_change_exit_x = zeros(nbins,1); % average change in exit rate by deciles of x
share_change_exit_x = zeros(nbins,1); % share of exit change by deciles of x

edges = linspace(0,1,nbins+1);


% x_c = 1
x_c = 1;d_c = 1;
if cum_x(x_c)>edges(d_c+1)
    % distrib pdf_x(x_c) to bins d_c and d_c+1 such that each bin gets
    % exactly the same measure of firms.

    % wgt is the weight on the left point (d_c), and 1-wgt is the weight on
    % the right point (d_c+1)
    wgt = edges(d_c+1)/cum_x(x_c);
else
    wgt = 1;
end

emp_decile_x(d_c) = emp_decile_x(d_c) + wgt*exitSS.emp_share_x(x_c);
output_decile_x(d_c) = output_decile_x(d_c) + wgt*exitSS.output_share_x(x_c);
ave_change_exit_x(d_c) = ave_change_exit_x(d_c) + ...
    wgt*(exit_grant_baseline.exit_rate_x(x_c) - exit_nogrant.exit_rate_x(x_c))*...
    (pdf_x(x_c))*nbins;
share_change_exit_x(d_c) = share_change_exit_x(d_c) + ...
    wgt*(exit_grant_baseline.exits_x(x_c) - exit_nogrant.exits_x(x_c));
if wgt<1
    emp_decile_x(d_c+1) = emp_decile_x(d_c+1) + (1-wgt)*exitSS.emp_share_x(x_c);
    output_decile_x(d_c+1) = output_decile_x(d_c+1) + (1-wgt)*exitSS.output_share_x(x_c);
    ave_change_exit_x(d_c+1) = ave_change_exit_x(d_c+1) + ...
        (1-wgt)*(exit_grant_baseline.exit_rate_x(x_c) - exit_nogrant.exit_rate_x(x_c))*...
        (pdf_x(x_c))*nbins;
    share_change_exit_x(d_c+1) = share_change_exit_x(d_c+1) + ...
        (1-wgt)*(exit_grant_baseline.exits_x(x_c) - exit_nogrant.exits_x(x_c));
    % update d_c
    d_c = d_c+1;
end
% x_c from 2 to nx
for x_c = 2:par.nx
    if cum_x(x_c)>edges(d_c+1) && cum_x(x_c-1)<edges(d_c+1)
        % distrib pdf_x(x_c) to bins d_c and d_c+1 such that each bin gets
        % exactly the same measure of firms.

        % wgt is the weight on the left point (d_c), and 1-wgt is the weight on
        % the right point (d_c+1)
        wgt = (edges(d_c+1)-cum_x(x_c-1))/(cum_x(x_c)-cum_x(x_c-1));
    else
        wgt = 1;
    end
    

    emp_decile_x(d_c) = emp_decile_x(d_c) + wgt*exitSS.emp_share_x(x_c);
    output_decile_x(d_c) = output_decile_x(d_c) + wgt*exitSS.output_share_x(x_c);
    ave_change_exit_x(d_c) = ave_change_exit_x(d_c) + ...
        wgt*(exit_grant_baseline.exit_rate_x(x_c) - exit_nogrant.exit_rate_x(x_c))*...
        (pdf_x(x_c))*nbins;
    share_change_exit_x(d_c) = share_change_exit_x(d_c) + ...
        wgt*(exit_grant_baseline.exits_x(x_c) - exit_nogrant.exits_x(x_c));

   
    if wgt<1
        emp_decile_x(d_c+1) = emp_decile_x(d_c+1) + (1-wgt)*exitSS.emp_share_x(x_c);
        output_decile_x(d_c+1) = output_decile_x(d_c+1) + (1-wgt)*exitSS.output_share_x(x_c);
        ave_change_exit_x(d_c+1) = ave_change_exit_x(d_c+1) + ...
            (1-wgt)*(exit_grant_baseline.exit_rate_x(x_c) - exit_nogrant.exit_rate_x(x_c))*...
            (pdf_x(x_c))*nbins;
        share_change_exit_x(d_c+1) = share_change_exit_x(d_c+1) + ...
            (1-wgt)*(exit_grant_baseline.exits_x(x_c) - exit_nogrant.exits_x(x_c));
         % update d_c
        d_c = d_c+1;
    end

end
share_change_exit_x = share_change_exit_x/sum(share_change_exit_x);


end % end function <fun_deciles_x>