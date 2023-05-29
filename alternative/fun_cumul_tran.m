function [sr_aveimp,cum_aveimp,ratio_aveimp] = fun_cumul_tran(irf,T_sr)
% Input: 
%   irf: struct with impulse responses in % relative to the steady state
%   T_sr: Short run length in quarters
% Outputs:
%   sr_aveimp: short-run average impulse response (in % relative to the steady state)
%   cum_aveimp: cumulative average impulse response (in % relative to the steady state)


% Cumulative impulse response under different policy environment 
%   and short- vs long-run cumulative responses.

% T_tran: length of transition
T_tran = length(irf.C_agg);


varNames = {'C_agg';'KL_ratio';'q';'w';'rental';
            'K_agg';'K_all';'K_small_owned';'K_corp';
            'mass_small';
            'Y_agg';'Y_corp';'output_small';
            'L_agg';'L_small';'L_corp';
            'entry'; 'exit';
            'InvK';'InvK_all'};

sr_aveimp      = struct();
cum_aveimp     = struct();
ratio_aveimp   = struct();
for ii = 1:numel(varNames)
    varii = varNames{ii};
    sr_aveimp.(varii)    = sum(irf.(varii)(1:T_sr))/T_sr;
    cum_aveimp.(varii)   = sum(irf.(varii))/T_tran;
    ratio_aveimp.(varii) = sum(irf.(varii)(1:T_sr))/sum(irf.(varii));
end

end % function fun_cumul_tran

