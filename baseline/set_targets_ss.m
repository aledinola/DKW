function [targetNames,targetNames_long,calibWeights,data_mom] = set_targets_ss()

%% Load data moments for the steady state
% Source: 'data_moments\data_moments.txt'
targetNames = ...
    {'avefirmsize';
    'avefirmsize_age0';
    'empshare_small';
    %'revshare_small';
    'exitrate';
    'fixedcost_to_rev';
    'autocorr_emp';
    'ave_work';
    'exitrate_0_9';
    'scor_invrate';	
    'hasNetDebt';
    'freq_lumpinv';
    'frac_exit_forced';
    'frac_exit_vol';
    'firmshare_0_9';
    'firmshare_10_19';
    'firmshare_20_99';
    'firmshare_100_499';
    'empshare_0_9';
    'empshare_10_19';
    'empshare_20_99';
    'empshare_100_499'}; 
% Moments not added to the above list:
% autocorr_logrev,fixedcost_to_payroll

targetNames_long = {'Average employment in small firms';
    'Average employment, age 0';
    'Small firm share of employment';
    'Small firm exit rate';
    'Fixed expense to revenue ratio';
    'Autocorr. employment';
    'Time spent in market work';
    'Exit rate, emp. size 0 to 9';
    'serial corr. investment rate';
    'Share of firms with debt';
    'Freq. positive lumpy investment';
    'Share forced exit';
    'Share voluntary exit';
    'firmshare0to9';
    'firmshare10to19';
    'firmshare20to99';
    'firmshare100to499';
    'empshare0to9';
    'empshare10to19';
    'empshare20to99';
    'empshare100to499'};

if ~isequal(numel(targetNames),numel(targetNames_long))
    error("character arrays <targetNames> and <targetNames_long> MUST have the same number of elements")
end

%% Set weights for calibration
calibWeights.avefirmsize    = 300;
calibWeights.empshare_small = 5;
calibWeights.exitrate       = 150;
calibWeights.frac_exit_forced = 0;
calibWeights.frac_exit_vol = 0;
calibWeights.avefirmsize_age0 = 300;
calibWeights.fixedcost_to_rev     = 5;
calibWeights.autocorr_emp   = 1000;
calibWeights.ave_work   = 100;
calibWeights.hasNetDebt = 0;
calibWeights.firmshare_0_9     = 1;
calibWeights.firmshare_10_19   = 1;
calibWeights.firmshare_20_99   = 1;
calibWeights.firmshare_100_499 = 1;
calibWeights.empshare_0_9     = 2;
calibWeights.empshare_10_19   = 2;
calibWeights.empshare_20_99   = 1;
calibWeights.empshare_100_499 = 1;
calibWeights.exitrate_0_9     = 30;
calibWeights.scor_invrate = 50;
calibWeights.freq_lumpinv = 5;

filename = fullfile('..','data_moments','data_moments.txt');
data_mom = read_data_targets(filename); %data_mom is a structure

end %end function "set_targets_ss"

