function [] = txt_export_tran(par,model_mom_trans,agg_tran,TabDir)

% MANUAL CALIBRATION OF THE TRANSITION SHOCKS
% txt_export_tran writes on txt file parameter values for shocks and
% targets


% Write the same in txt file
fid=fopen(fullfile(TabDir,'targets_transition_manual.txt'),'wt');  % overwrite

fprintf(fid,'PARAMETER \n');
fprintf(fid,'eta_i         :    %f \n', par.eta_i);
fprintf(fid,'v_small       :    %f \n', par.v_small);
fprintf(fid,'v_small_unimp :    %f \n', par.v_small_unimp);
fprintf(fid,'v_corp        :    %f \n', par.v_corp);
fprintf(fid,'util_shift    :    %f \n', par.util_shift);
fprintf(fid,'lsupply_shift :    %f \n', par.lsupply_shift);
fprintf(fid,'lambda_shift  :    %f \n', par.lambda_shift);
fprintf(fid,'mass_shift  :    %f \n', par.mass_shift);
fprintf(fid,'rho_shock     :    %f \n', par.rho_shock);

fprintf(fid,' \n');

fprintf(fid,'MOMENT \n');
fprintf(fid,'Change in GDP q1                   :    %f \n', model_mom_trans(1,1));
fprintf(fid,'Change in GDP q2                   :    %f \n', model_mom_trans(1,2));
fprintf(fid,'Change in consumption              :    %f \n', model_mom_trans(2,1));
fprintf(fid,'Change in investment               :    %f \n', model_mom_trans(3,1));
fprintf(fid,'Change in small firms output       :    %f \n', model_mom_trans(4,1));
fprintf(fid,'Change in employment q1            :    %f \n', model_mom_trans(5,1));
fprintf(fid,'Change in employment q2            :    %f \n', model_mom_trans(5,2));
fprintf(fid,'Change in employment small q1      :    %f \n', model_mom_trans(6,1));
fprintf(fid,'Change in employment small q2      :    %f \n', model_mom_trans(6,2));
fprintf(fid,'Change in employment corp q1       :    %f \n', model_mom_trans(7,1));
fprintf(fid,'Change in employment corp q2       :    %f \n', model_mom_trans(7,2));
fprintf(fid,'Change in small firm exit q1       :    %f \n', model_mom_trans(8,1));
fprintf(fid,'Change in small firm exit, annual  :    %f \n', model_mom_trans(9,1));
fprintf(fid,'Change in small firm entry, q1  :    %f \n', model_mom_trans(10,1));
fprintf(fid,'Total grant amount                 :    %f \n', agg_tran.tot_grant);

fclose(fid);




end % end function