function [] = append_results_txt(obj_smm,param_struct,calib_names,data_mom_struct,model_mom_struct,calib_weights_struct,mom_names,TxtDir)

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

% Convert structure with parameters into a vector
param_vec = zeros(numel(calib_names),1);
for i = 1:numel(calib_names) 
    param_vec(i) = param_struct.(calib_names{i});
end

% Convert data_moments into vector
data_mom_vec = zeros(numel(mom_names),1); % data moments
for i = 1:numel(mom_names)
    data_mom_vec(i) = data_mom_struct.(mom_names{i});
end

% Convert model_moments into vector
model_mom_vec = zeros(numel(mom_names),1); % model moments
for i = 1:numel(mom_names)
    model_mom_vec(i) = model_mom_struct.(mom_names{i});
end

% Convert calibration weigths into vector
calibWeights_vec = zeros(numel(mom_names),1); % calibration weights
for i = 1:numel(mom_names)
    calibWeights_vec(i) = calib_weights_struct.(mom_names{i});
end

dev = zeros(numel(mom_names),1); % abs deviation
for i = 1:numel(mom_names)
    dev(i) = calibWeights_vec(i)*((data_mom_vec(i)-model_mom_vec(i))/data_mom_vec(i))^2;
end

%% Open file
FID = fopen(fullfile(TxtDir,'results_sofar.txt'),'a+');

%% Append parameters to txt file
width = max(cellfun('length', calib_names));
fprintf(FID,'--------------------------------------------------  \n');
fprintf(FID,' Parameter  Value  \n');
fprintf(FID,'--------------------------------------------------  \n');
for i = 1:numel(calib_names)
    fprintf(FID,'%-*s   %-8.16f \n',width,calib_names{i},param_vec(i));
end
fprintf(' \n')

%% Append model fit results to txt file
width = max(cellfun('length', mom_names)); 
fprintf(FID,'--------------------------------------------------  \n');
fprintf(FID,'%-*s %-s         %-s \n',width,'Moment','Data','Model');
fprintf(FID,'--------------------------------------------------  \n');
for i = 1:numel(mom_names)
    fprintf(FID,'%-*s    %-8.4f  %-8.4f  \n',width,mom_names{i},data_mom_vec(i),model_mom_vec(i));
end
fprintf(FID,'--------------------------------------------------  \n');
fprintf(FID,'obj_smm = %-8.4f \n',obj_smm);
fprintf(FID,'=======================================================  \n');

%% Close file
fclose(FID);

end %END FUNCTION

