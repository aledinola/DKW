function [data_mom] = read_data_targets(filename)
%Purpose: read data targets (names and values) from txt file filename

data_mom = struct();

FID = fopen(filename);
C = textscan(FID,'%s %f');
fclose(FID);

names  = C{1};
values = C{2};

for i=1:numel(names)
   data_mom.(names{i}) = values(i); 
end

% Manual conversion
data_mom.firmshare_0_9 = data_mom.firmshare_0_4+data_mom.firmshare_5_9;
data_mom.empshare_0_9 = data_mom.empshare_0_4 + data_mom.empshare_5_9;
data_mom.empshare_small = data_mom.empshare_small_imp+data_mom.empshare_small_unimp;
data_mom.revshare_small = data_mom.revshare_small_imp+data_mom.revshare_small_unimp;
data_mom.ave_work       = 0.33;
data_mom.frac_exit_forced = 1.0;
data_mom.frac_exit_vol = 1.0;

data_mom.debt_asset_all_entrants = data_mom.debt_asset_all/data_mom.debt_asset_entrants;


% Convert annual rates to quarterly rates
% Transition rates in the data are based on annual data:
% - exit rate, job creation and destruction rates are from BDS
% - autocorrepation of emp from KFS

data_mom.exitrate =  1- (1-data_mom.exitrate)^(1/4);
data_mom.exitrate_0_9 =  1- (1-data_mom.exitrate_0_9)^(1/4);
data_mom.exitrate_10_19 =  1- (1-data_mom.exitrate_10_19)^(1/4);
data_mom.exitrate_20_99 =  1- (1-data_mom.exitrate_20_99)^(1/4);
data_mom.exitrate_100_499 =  1- (1-data_mom.exitrate_100_499)^(1/4);
data_mom.jcr = data_mom.jcr/4;
data_mom.jdr = data_mom.jdr/4;
data_mom.autocorr_emp = data_mom.autocorr_emp^(1/4);

end %end function

