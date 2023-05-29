function [] = append_tran_txt(distance,x_in,data_mom_trans,model_mom_trans,calibWeightsTran,InpDir)

%{
DESCRIPTION:
Write intermediate results (parameter values, model fit and distance b/w
model moments and data moments) to a txt file. Especially useful when
estimating parameters with a global optimization routine.

INPUTS:
- param_struct: structure with all parameters
- calib_names: cell array of characters with names of estimated parameters 
- mom_names: cell array of characters with names of targeted moments
- data_mom_struct:  structure with data moments
- model_mom_struct: structure with model moments
- calib_weights_struct: structure with calibration weigths
- TxtDir:    Folder where you want to append txt file
%}

%% Open file
FID = fopen(fullfile(InpDir,'results_tran_sofar.txt'),'a+');

%% Append parameters to txt file
fprintf(FID,"==================================== \n");
fprintf(FID,"eta_i:       %8.6f  \n",x_in(1));
fprintf(FID,"v_corp:        %8.6f  \n",x_in(2));
fprintf(FID,"util_shift:    %8.6f  \n",x_in(3));
fprintf(FID,"lsupply_shift: %8.6f  \n",x_in(4));
fprintf(FID,"rho_shock:     %8.6f  \n",x_in(5));
fprintf(FID," \n");

%% Append transition fit to txt file
width = length('Drop in small firms exit rate:')+3;
fprintf(FID,"  \n");
fprintf(FID,"%-*s  %-s     %-s     %-s     \n",width,"Description","Data","Model","Weight");
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in GDP q1:",       data_mom_trans(1,1),model_mom_trans(1,1),calibWeightsTran(1,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in GDP q2:",       data_mom_trans(1,2),model_mom_trans(1,2),calibWeightsTran(1,2));
%fprintf("%-*s %-8.4f  %-8.4f  \n", width,"Drop in GDP q3:",               data_mom_trans(1,3),model_mom_trans(1,3))
%fprintf("%-*s %-8.4f  %-8.4f  \n", width,"Drop in GDP q4:",               data_mom_trans(1,4),model_mom_trans(1,4))
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f\n", width,"Drop in consumption:",          data_mom_trans(2,1),model_mom_trans(2,1),calibWeightsTran(2,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in investment:",         data_mom_trans(3,1),model_mom_trans(3,1),calibWeightsTran(3,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in small firms output:", data_mom_trans(4,1),model_mom_trans(4,1),calibWeightsTran(4,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment q1:",      data_mom_trans(5,1),model_mom_trans(5,1),calibWeightsTran(5,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment q2:",      data_mom_trans(5,2),model_mom_trans(5,2),calibWeightsTran(5,2));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment small q1:",data_mom_trans(6,1),model_mom_trans(6,1),calibWeightsTran(6,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment small q2:",data_mom_trans(6,2),model_mom_trans(6,2),calibWeightsTran(6,2));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment corp q1:", data_mom_trans(7,1),model_mom_trans(7,1),calibWeightsTran(7,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in Employment corp q2:", data_mom_trans(7,2),model_mom_trans(7,2),calibWeightsTran(7,2));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Drop in small firms exit rate:", data_mom_trans(8,1),model_mom_trans(8,1),calibWeightsTran(8,1));
fprintf(FID,"%-*s  %-8.4f  %-8.4f  %-8.4f  \n", width,"Change in annual exit rate:", data_mom_trans(9,1),model_mom_trans(9,1),calibWeightsTran(9,1));

fprintf(FID,"%-*s %-8.6f          \n",width,"Distance:",distance);

fprintf(FID,"  \n");

%% Close file
fclose(FID);

%% Old stuff
% % Convert structure with parameters into a vector
% param_vec = zeros(numel(calib_names),1);
% for i = 1:numel(calib_names) 
%     param_vec(i) = param_struct.(calib_names{i});
% end
% 
% % Convert data_moments into vector
% data_mom_vec = zeros(numel(mom_names),1); % data moments
% for i = 1:numel(mom_names)
%     data_mom_vec(i) = data_mom_struct.(mom_names{i});
% end
% 
% % Convert model_moments into vector
% model_mom_vec = zeros(numel(mom_names),1); % model moments
% for i = 1:numel(mom_names)
%     model_mom_vec(i) = model_mom_struct.(mom_names{i});
% end
% 
% % Convert calibration weigths into vector
% calibWeights_vec = zeros(numel(mom_names),1); % calibration weights
% for i = 1:numel(mom_names)
%     calibWeights_vec(i) = calib_weights_struct.(mom_names{i});
% end
% 
% dev = zeros(numel(mom_names),1); % abs deviation
% for i = 1:numel(mom_names)
%     dev(i) = calibWeights_vec(i)*((data_mom_vec(i)-model_mom_vec(i))/data_mom_vec(i))^2;
% end
% 
% %% Open file
% FID = fopen(fullfile(TxtDir,'results_sofar.txt'),'a+');
% 
% %% Append parameters to txt file
% width = max(cellfun('length', calib_names));
% fprintf(FID,'--------------------------------------------------  \n');
% fprintf(FID,' Parameter  Value  \n');
% fprintf(FID,'--------------------------------------------------  \n');
% for i = 1:numel(calib_names)
%     fprintf(FID,'%-*s   %-8.16f \n',width,calib_names{i},param_vec(i));
% end
% fprintf(' \n')
% 
% %% Append model fit results to txt file
% width = max(cellfun('length', mom_names)); 
% fprintf(FID,'--------------------------------------------------  \n');
% fprintf(FID,'%-*s %-s         %-s \n',width,'Moment','Data','Model');
% fprintf(FID,'--------------------------------------------------  \n');
% for i = 1:numel(mom_names)
%     fprintf(FID,'%-*s    %-8.4f  %-8.4f  \n',width,mom_names{i},data_mom_vec(i),model_mom_vec(i));
% end
% fprintf(FID,'--------------------------------------------------  \n');
% fprintf(FID,'obj_smm = %-8.4f \n',obj_smm);
% fprintf(FID,'=======================================================  \n');
% 
% %% Close file
% fclose(FID);

end %END FUNCTION

