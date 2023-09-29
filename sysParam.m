%% FUNCTION - Store System Parameters
function param = sysParam()

    % requirement
    param.vrq = 10218;
    param.mPL = 24; % kg
    param.Mearth = 5.9722e+24;
    param.Rearth = 6378100;
    param.release_vel = 250.786;
    param.release_alt = 12192;
    param.release_temp = 216.65;
    param.orb_alt = 251460;
    param.orb_vel = 7754.1;
    param.bounding_box_side = 24*0.0254;
    param.bounding_box_leng = 240*0.0254;
    param.bounding_box_volu = param.bounding_box_side^2*param.bounding_box_leng;
    param.bounding_mass = 2265;
    param.scram_ceiling = 50000;
    param.scram_ceiling_temp = 270.65;
    param.scram_mach = 8;
    param.scram_solid_boost_mach = 3;

    % universal parameters
    param.g0 = 9.80655;

    param.is_scram = true;
    param.is_scram_solid_boost = false;
    % first stage
    param.Isp_stg1 = 1200;
    param.sigma_stg1 = 0.4;
    param.density_stg1 = 800;
    
    % second stage
    param.Isp_stg2 = 270;
    param.sigma_stg2 = 0.1;
    param.density_stg2 = 1370;

    % third stage
    param.Isp_stg3 = 320;
    param.sigma_stg3 = 0.1;
    param.density_stg3 = 2040;

    % dv losses
    param.dragloss = 65;
    param.proploss = 150;
    param.steeloss = 360;
    param.manvloss = 50;
    