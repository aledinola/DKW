%% This script generates latex tables

if ~(isfolder(TabDir))
    mkdir(TabDir)
end

fid = fopen(fullfile(TabDir, 'tableI.tex'), 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Steady State Results} \n');
%fprintf(fid,'\\label{tab:results_us} \n');
fprintf(fid, '\\begin{tabular}{lc} \\hline \n');
fprintf(fid, '\\hline \n');
%fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,'Moment                     & Value  \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid,"Measure of active firms:   &  %8.3f     \\\\ \n",Mactive );
fprintf(fid,"Measure of entrants:       &  %8.3f     \\\\ \n",Mentr );
fprintf(fid,"Aggregate premises:        &  %8.3f     \\\\ \n",P_agg );
fprintf(fid,"Investment into premises:  &  %8.3f     \\\\ \n",InvP  );
fprintf(fid,"Aggregate consumption:     &  %8.3f     \\\\ \n",C_agg );
fprintf(fid,"Capital (corporate):       &  %8.3f     \\\\ \n",K_corp );
fprintf(fid,"Employment (corporate):    &  %8.3f     \\\\ \n",L_corp );
fprintf(fid,"Employment (total):        &  %8.3f     \\\\ \n",N_agg );
fprintf(fid,"Output (corporate):        &  %8.3f     \\\\ \n",Y_corp );
fprintf(fid,"Output (small firms):      &  %8.3f     \\\\ \n",Y_small );
fprintf(fid,"Output (total):            &  %8.3f     \\\\ \n",Y_small+Y_corp );

fprintf(fid,'\\hline \n \\end{tabular} \n');
% End tabular environment

fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);
