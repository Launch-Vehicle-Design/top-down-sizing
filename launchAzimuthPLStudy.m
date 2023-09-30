clear; clc; close all

load("OptimalSolution.mat")

%% defined system
% general parameters
op_site_lat = 0:90; op_orb_inc = 0:90;
[OP_SITE_LAT,OP_ORB_INC] = meshgrid(op_site_lat,op_orb_inc);
cos_ratio = cos(OP_ORB_INC/180*pi)./cos(OP_SITE_LAT/180*pi);
cos_ratio(cos_ratio > 1) = 1;
OP_INER_AZ_ASD = asin(cos_ratio); OP_INER_AZ_DSD = pi-asin(cos_ratio);
OP_INER_AZ_ASD(OP_SITE_LAT>OP_ORB_INC) = nan;
OP_INER_AZ_DSD(OP_SITE_LAT>OP_ORB_INC) = nan;
OP_SITE_VEL = 465.1*cos(OP_SITE_LAT/180*pi);

cal_losses = 2272.52; cal_reserved = 200.54;
g0 = param.g0; Isp1 = param.Isp_stg1;
Isp2 = param.Isp_stg2; Isp3 = param.Isp_stg3;

% unmatched site orbit launch dv budget
REQ_ROT_VEL = sqrt((param.orb_vel*sin(OP_INER_AZ_ASD)-OP_SITE_VEL).^2+(param.orb_vel*cos(OP_INER_AZ_ASD)).^2);
REQ_UNMATCH_DV_BUDGET = REQ_ROT_VEL+cal_losses+cal_reserved-param.release_vel;

% 1st stage
m01 = optimal_3stg(ind(1,"m0"));
mp1 = optimal_3stg(ind(1,"mp"));
ms1 = optimal_3stg(ind(1,"ms"));
m01_unw = m01-param.mPL;

% 2nd stage
m02 = optimal_3stg(ind(2,"m0"));
mp2 = optimal_3stg(ind(2,"mp"));
ms2 = optimal_3stg(ind(2,"ms"));
m02_unw = m02-param.mPL;

% 3rd stage
m03 = optimal_3stg(ind(3,"m0"));
mp3 = optimal_3stg(ind(3,"mp"));
ms3 = optimal_3stg(ind(3,"ms"));
m03_unw = m03-param.mPL;

% 3 stg input check - if total dv is delivered for std PL
deliver_dv = g0*(Isp1*log(m01/(m01-mp1))+Isp2*log(m02/(m02-mp2))+Isp3*log(m03/(m03-mp3)));
if deliver_dv < param.vrq
    disp("ERROR - System with Insufficient dv"); return
end

% iteratively solve for unmatch cases
dvf = @(mPL) g0*(Isp1*log((m01_unw+mPL)/(m01_unw+mPL-mp1))+ ...
    Isp2*log((m02_unw+mPL)/(m02_unw+mPL-mp2))+ ...
    Isp3*log((m03_unw+mPL)/(m03_unw+mPL-mp3)));

MPL_SOLVE = nan(size(OP_SITE_LAT));
for lat = 1:size(OP_SITE_LAT,2)
    for inc = size(OP_SITE_LAT,1):-1:lat
        % REQ_UNMATCH_DV_BUDGET(inc,lat)
        ddv = @(mPL) dvf(mPL)-REQ_UNMATCH_DV_BUDGET(inc,lat);
        MPL_SOLVE(inc,lat) = fzero(ddv,param.mPL);
    end
end

figure; surf(OP_SITE_LAT,OP_ORB_INC,MPL_SOLVE);
xlabel("Launch site latitude"); ylabel("Desired orbit inclination");
zlabel("Payload capacity");

%% Helper function ind for 3 stages
function i = ind(stage,mass_type)
    ind_m0 = 1; ind_ms = 2; ind_mp = 3; num_param = 5;
    if mass_type == "m0"
        i = (3-stage)*num_param + ind_m0;
    elseif mass_type == "ms"
        i = (3-stage)*num_param + ind_ms;
    elseif mass_type == "mp"
        i = (3-stage)*num_param + ind_mp;
    else
        disp("ERROR - Not Recognizable Mass Type")
    end
end