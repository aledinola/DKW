function [] = make_table_ss(prices,agg,b_grid,do_tex,tabDir,filename)

%{
DESCRIPTION:
Write steady-state results on the screen and (if do_tex=1) to a latex file

INPUTS:
- prices,agg,b_grid: Objects to write
- do_tex:            Flag 0/1 to write latex file
- tabDir:            Folder where you want to save the latex file
- filename:          Name of the latex file, must have .tex suffix

%}

if ~ischar(tabDir)
    error('Input tabDir in make_table_ss must be a character var')
end
if ~ischar(filename)
    error('Input filename in make_table_ss must be a character var')
end

if do_tex~=0 && do_tex~=1
    error('do_tex is out of range: must be either 0 or 1')
end

nk = size(b_grid,1);

fprintf("  \n")
disp("--------------------------------------------")
disp("STEADY-STATE RESULTS")
disp("--------------------------------------------")

fprintf("Firms Discount factor q:   %10.16f  \n",prices.q )
fprintf("Wage rate:                 %10.16f  \n",prices.wage )
fprintf("Capital-labor ratio corp:  %10.16f  \n",prices.KL_ratio )

fprintf("Measure of active firms:   %10.16f  \n",agg.Mactive )
fprintf("Measure of entrants:       %10.16f  \n",agg.Mentr )
fprintf("Aggregate consumption:     %10.16f  \n",agg.C_agg )
fprintf("Capital (corporate):       %10.16f  \n",agg.K_corp )
fprintf("Employment (corporate):    %10.16f  \n",agg.L_corp )
fprintf("Employment (small firms):  %10.16f  \n",agg.L_small )
fprintf("Aggregate Capital:         %10.16f  \n",agg.K_agg )
fprintf("Employment (total):        %10.16f  \n",agg.L_agg )
fprintf("Output (corporate):        %10.16f  \n",agg.Y_corp )
fprintf("Output (small firms):      %10.16f  \n",agg.output_small )
fprintf("Total output:              %10.16f  \n",agg.Y_agg )
fprintf("liq:                       %10.16f  \n",agg.liq )
fprintf("entry:                     %10.16f  \n",agg.entry )
fprintf("Aggregate Investment:      %10.16f  \n",agg.InvK )
fprintf("btilde(1):                 %10.16f  \n",b_grid(1,1) ) %b_tilde has dim: (nk,1)
fprintf("btilde(nk):                %10.16f  \n",b_grid(nk,1) ) %b_tilde has dim: (nk,1)


fprintf("  \n")

% If desired, write the same results on a latex file
if do_tex==1
    
    FID = fopen(fullfile(tabDir,filename),'w');
    
    fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
    fprintf(FID,' \\hline \n');
    
    fprintf(FID,"Measure of active firms:   &  %8.4f  \\\\ \n",agg.Mactive );
    fprintf(FID,"Measure of entrants:       &  %8.4f  \\\\ \n",agg.Mentr );
    fprintf(FID,"Aggregate consumption:     &  %8.4f  \\\\ \n",agg.C_agg );
    fprintf(FID,"Capital (corporate):       &  %8.4f  \\\\ \n",agg.K_corp );
    fprintf(FID,"Employment (corporate):    &  %8.4f  \\\\ \n",agg.L_corp );
    fprintf(FID,"Aggregate Capital:         &  %8.4f  \\\\ \n",agg.K_agg );
    fprintf(FID,"Employment (total):        &  %8.4f  \\\\ \n",agg.L_agg );
    fprintf(FID,"Output (corporate):        &  %8.4f  \\\\ \n",agg.Y_corp );
    fprintf(FID,"Output (small firms):      &  %8.4f  \\\\ \n",agg.output_small );
    fprintf(FID,"Total output:              &  %8.4f  \\\\ \n",agg.Y_agg );
    fprintf(FID,"liq:                       &  %8.4f  \\\\ \n",agg.liq );
    fprintf(FID,"entry:                     &  %8.4f  \\\\ \n",agg.entry );
    fprintf(FID,"Aggregate Investment:      &  %8.4f  \\\\ \n",agg.InvK );
    fprintf(FID,"btilde(1):                 &  %8.4f  \\\\ \n",b_grid(1,1) ); %b_tilde has dim: (nk,1)
    fprintf(FID,"btilde(nk):                &  %8.4f  \\\\ \n",b_grid(nk,1) ); %b_tilde has dim: (nk,1)
    
    fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
    
    fclose(FID);
end

end % end function <make_table_ss>