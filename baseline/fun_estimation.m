function [obj_smm] = fun_estimation(model_mom,data_mom,targetNames,calibWeights)
% Inputs
% model_mom :: struct with model moments
% data_mom  :: struct with data moments
% targetNames :: character array with names of moments

% Convert model_mom struct into a vector
model_mom_vec = struct2vec(model_mom,targetNames);
data_mom_vec  = struct2vec(data_mom,targetNames);
calibWeights_vec = struct2vec(calibWeights,targetNames);

if length(model_mom_vec)~=length(data_mom_vec)
    error("Model and data targets are not compatible!")
end
if length(calibWeights_vec)~=length(data_mom_vec)
    error("Model and data targets are not compatible!")
end


obj_smm = 0;
for i=1:length(model_mom_vec)
    dist2 = ((model_mom_vec(i)-data_mom_vec(i))/data_mom_vec(i))^2;
    obj_smm = obj_smm+dist2*calibWeights_vec(i);
end

check = model_mom.empshare_small>1;
penalty = 1000*max(model_mom.empshare_small-1,0)^2;
if check
    obj_smm = obj_smm+penalty;
end

end %end function 

