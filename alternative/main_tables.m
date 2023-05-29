%% This script compiles results (tables) for the paper
% This script requires the following mat files:
%   - ss.mat, nogrant.mat, grant_baseline.mat, grant_targslim.mat 
% This script creates tables used in the writeup.
% loading saved results in mat files
% Tables created by this file:
%   - comp_parameters.tex: Table with Computational parameters
%   - exo_parameters.tex: Table with exogenous parameters
%   - steady_state.tex: Table with some steady-state values
%   - moments.tex: Table with model fit
%   - parameters.tex: Table with estimated parameters
%   - tran_shocks.tex: Table with shocks
%   - transition_moments.tex: table with transition moments
%   - cost_grants.tex: table with cost per job saved

clear;clc;close all

matNames = {'ss','nogrant','grant_baseline','grant_targslim'};
for ii = 1:numel(matNames)
    if ~isfile(fullfile('mat',[matNames{ii},'.mat']))
        warning('MAT file ''%s'' is missing \n',matNames{ii})
    end
end

%% Make tables for steady-state

% Load steady-state results
if isfile(fullfile('mat','ss.mat'))
    load(fullfile('mat','ss.mat'))
else
    error('File ss.mat does not exist!')
end

% Table with estimated parameters (the same for both economies)
% mystruct2table(pstruct,pnames,description,dispNames,header,tex,tabDir,filename)
mystruct2table(par,calibNames,description,dispNames,{'Parameter','Value'},1,fullfile('tables'),'parameters.tex');

% Table with model fit
mystruct2table_mom(data_mom,model_mom,targetNames,calibWeights,targetNames_long,1,fullfile('tables'),'moments.tex');
