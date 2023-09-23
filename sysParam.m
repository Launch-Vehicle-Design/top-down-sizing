%% FUNCTION - Store System Parameters
function param = sysParam()

    % requirement
    param.vrq = 10218;
    param.mPL = 24; % kg
    param.Mearth = 5.9722e+24;
    param.Rearth = 6378100;
    param.release_vel = 250.786;
    param.release_alt = 12192;
    param.orb_vel = 251460;
    param.orb_alt = 7754.1;

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

    % dv losses
    param.dragloss = 65;
    param.proploss = 150;
    param.steeloss = 360;
    param.manvloss = 50;
    