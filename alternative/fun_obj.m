function [obj_smm,sol,agg,b_grid,distribS,prices,model_mom,flag_ss,par] = fun_obj(guess,par,bounds,calibNames,data_mom,targetNames,calibWeights,description,dispNames,targetNames_long)
%{
This function takes a vector of parameter values "guess" as input, solves the
model and generates the distance b/w model moments and data moments
INPUTS:
- guess:      vector with parameters to be calibrated
- par:        structure with parameters
- bounds:     structure with bounds for calibration parameters. Each field
              is 1*2 vector
- calibNames: cell array of strings with names of parameters to be
              calibrated/estimated
- data_mom:         vector with data moments
- targetNames:      cell array of strings
- calibWeights:     numerical array
- description:      cell array of strings
- dispNames:        cell array of strings (latex names of parameters)
- targetNames_long: cell array of strings

OUTPUTS:
- obj_smm:   distance between model targets and data targets
%}

global obj_smm_best %#ok<GVMIS> 

% Input checks:
if ~isstruct(par)
    error("input par in <fun_obj> must be a struct!")
end
if ~isstruct(bounds)
    error("input bounds in <fun_obj> must be a struct!")
end
if isrow(guess)
   error("guess must be a column vector!") 
end
if size(guess,1)~=size(calibNames,1)
    error("guess and calibNames must have the same no of elements!")
end
if  length(fieldnames(bounds))~=size(guess,1)
    error("guess and bounds do not match!")
end
if ~iscell(targetNames)
    error("targetNames must be a cell array")
end
if ~iscell(description)
    error("description must be a cell array")
end
if ~iscell(targetNames_long)
    error("targetNames_long must be a cell array")
end
if ~iscell(dispNames)
    error("dispNames must be a cell array")
end

% Add parameters to be estimated to the structure PAR
[par] = vec2struct(guess,calibNames,par);

%Convert bounds from struct to vector
bounds_vec = bounds2vec(bounds,calibNames);

lbounds_vec = bounds_vec(:,1); %LOWER bounds
ubounds_vec = bounds_vec(:,2); %UPPER bounds

%check if parameters are within lower bounds
if any(guess<lbounds_vec)
    warning('Lower bounds on parameters violated!')
    obj_smm = 10000; %1e+10
    sol=[];agg=[];b_grid=[];distribS=[];prices=[];model_mom=[];flag_ss=[];par=[];
    return
end
%check if parameters are within upper bounds
if any(guess>ubounds_vec)
    warning('Upper bounds on parameters violated!')
    obj_smm = 10000;
    sol=[];agg=[];b_grid=[];distribS=[];prices=[];model_mom=[];flag_ss=[];par=[];
    return
end

%% Solve for the steady-state
if par.verbose>=2
    fprintf('Start model solution... \n')
    fprintf(' \n')
    fprintf(' \n')
end

%tic
% fun_steady_state solves the steady-state of the model and produces
% solution objects and model_mom
[sol,agg,b_grid,distribS,prices,model_mom,flag_ss,par] = fun_steady_state(par);
%toc

if flag_ss<0
    % Either VF or distrib did not converge
    warning('Some error occurred in fun_steady_state')
    obj_smm = 10000;
    sol=[];agg=[];b_grid=[];distribS=[];prices=[];model_mom=[];flag_ss=[];par=[];
    return
end

if agg.K_corp<0
    warning("Capital in corporate sector is negative!")
end

%% Table with estimated parameters
if par.do_table==1
    mystruct2table(par,calibNames,description,dispNames,{'Parameter','Value'},par.do_tex,par.TabDir,'parameters.tex');
end

%% Table with targeted moments: model vs data
if par.do_table==1
    mystruct2table_mom(data_mom,model_mom,targetNames,calibWeights,targetNames_long,par.do_tex,par.TabDir,'moments.tex');
end

%% Plots steady-state

if par.do_plots_ss==1
    make_plots_ss(sol,distribS,b_grid,par);
end

%% Assign distance for output
obj_smm = fun_estimation(model_mom,data_mom,targetNames,calibWeights);

if isnan(obj_smm)
    obj_smm = 10000;
end

fprintf('\n');
fprintf("obj_smm:      %f  \n",obj_smm )
fprintf('=============================================================\n');

if obj_smm<obj_smm_best
    obj_smm_best = obj_smm;
    % append results to txt file
    append_results_txt(obj_smm,par,calibNames,data_mom,model_mom,calibWeights,targetNames,par.InpDir);
end

end %end function

