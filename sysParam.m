%% FUNCTION - Store System Parameters
function param = sysParam()

    % requirement
    param.vrq = 10602;
    param.mPL = 24; % kg

    % universal parameters
    param.g0 = 9.80655;

    % first stage
    param.Isp_stg1 = 260;
    param.sigma_stg1 = 0.05;
    
    % second stage
    param.Isp_stg2 = 285.6;
    param.sigma_stg2 = 0.08;

    % third stage
    param.Isp_stg3 = 285.6;
    param.sigma_stg3 = 0.1;
    