function [targetNames,targetNames_long,calibWeights,data_mom] = set_targets_ss()
% FOLDER: 'alternative'

%% Load data moments for the steady state
% Source: shared forlder, data, moments, data_moments.txt
% Remove 'autocorr_emp'
targetNames = ...
    {'jcr';
    'jdr';
    'avefirmsize';
    'empshare_small';
    'exitrate';
    'avefirmsize_age0';
    'fixedcost_to_rev';
    'autocorr_emp';
    'ave_work';
    'exitrate_0_9';
    'scor_invrate'}; 

targetNames_long = {'Job creation rate';
    'Job destruction rate';  
    'Average employment in small firms';
    'Small firm share of employment';
    'Small firm exit rate';
    'Average employment, age 0';
    'Fixed expense to revenue ratio';
    'Autocorr. employment';
    'Time spent in market work';
    'Exit rate, emp. size 0-9';
    'Serial corr. investment rate'};

if ~isequal(numel(targetNames),numel(targetNames_long))
    error("character arrays <targetNames> and <targetNames_long> MUST have the same number of elements")
end

calibWeights.jcr = 0.2;
calibWeights.jdr = 0;
calibWeights.avefirmsize    = 300;
calibWeights.empshare_small = 5;
calibWeights.exitrate       = 150;
calibWeights.avefirmsize_age0 = 300;
calibWeights.fixedcost_to_rev     = 5;
calibWeights.autocorr_emp   = 1000;
calibWeights.ave_work   = 100;
calibWeights.exitrate_0_9     = 30;
calibWeights.scor_invrate = 50;

filename = fullfile('..','data_moments','data_moments.txt');
data_mom = read_data_targets(filename); %data_mom is a structure


end %end function "set_targets_ss"

