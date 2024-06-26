%% FUNCTION - Store System Parameters
function param = sysParam()

    % requirement
    param.mPL = 46; % kg
    param.mPL_lowerBound = 24;
    param.mPL_highBound = 45.4;
    param.G = 0.0000000000667428;
    param.Mearth = 5.9722e+24;
    param.Rearth = 6378100;
    param.release_lat = 80/180*pi;
    param.release_vel = 250.786;
    param.release_alt = 12192;
    param.release_temp = 216.65;
    param.orb_alt = 251460;
    param.orb_vel = 7754.1;
    param.bounding_box_side = 24*0.0254;
    param.bounding_box_leng = 240*0.0254;
    param.bounding_box_volu = param.bounding_box_side^2*param.bounding_box_leng;
    param.bounding_mass = 2265;
    param.scram_ceiling = 75000;
    param.scram_ceiling_temp = 270.65;
    param.scram_mach = 9;
    param.scram_solid_boost_mach = 5;

    % universal parameters
    param.g0 = 9.80665;

    param.is_scram = true;
    param.is_scram_solid_boost = false;
    % first stage
    param.Isp_stg1 = 290;
    param.sigma_stg1 = 0.06;
    param.density_stg1 = 2040;
    
    % second stage
    param.Isp_stg2 = 290;
    param.sigma_stg2 = 0.08;
    param.density_stg2 = 2040;

    % third stage
    param.Isp_stg3 = 290;
    param.sigma_stg3 = 0.1;
    param.density_stg3 = 2040;

    % dv losses
    param.gravloss = 1050;
    param.dragloss = 40;
    param.proploss = 150;
    param.steeloss = 200;
    reserv_perc = 0.02;
    param.dvreserv = (param.gravloss+param.dragloss+param.proploss+param.steeloss+param.orb_vel)*reserv_perc;

    % total dv
    param.vrq = param.dvreserv/reserv_perc*(1+reserv_perc)-param.release_vel+465.1*cos(param.release_lat);
    