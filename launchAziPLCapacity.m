function [MPL,op_site_lat,op_orb_inc] = launchAziPLCapacity(optimal, param, is_3stg, title_name)

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

plot_imperial = true;
if exist("title","var")
    title_name = "Impact of Launch Location and Inclination on PL Capacity";
end

%% defined system
% general parameters
op_site_lat = 0:0.2:90; op_orb_inc = 0:0.2:180;
[OP_SITE_LAT,OP_ORB_INC] = meshgrid(op_site_lat,op_orb_inc);
cos_ratio = cos(OP_ORB_INC/180*pi)./cos(OP_SITE_LAT/180*pi);
if sum(cos_ratio(cos_ratio > 1 & cos_ratio < 1+0.0001)) ~= 0
    cos_ratio(cos_ratio > 1) = 1;
end
if sum(cos_ratio(cos_ratio > 1)) ~= 0
    cos_ratio(cos_ratio > 1) = nan;
end
if sum(cos_ratio(cos_ratio < -1 & cos_ratio > -1-0.0001)) ~= 0
    cos_ratio(cos_ratio < -1 & cos_ratio > -1-0.0001) = -1;
end
if sum(cos_ratio(cos_ratio < -1)) ~= 0
    cos_ratio(cos_ratio < -1) = nan;
end
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

if is_3stg
    % 1st stage
    m01 = optimal(ind(3,1,"m0"));
    mp1 = optimal(ind(3,1,"mp"));
    ms1 = optimal(ind(3,1,"ms"));
    m01_unw = m01-param.mPL;
    
    % 2nd stage
    m02 = optimal(ind(3,2,"m0"));
    mp2 = optimal(ind(3,2,"mp"));
    ms2 = optimal(ind(3,2,"ms"));
    m02_unw = m02-param.mPL;
    % 3rd stage
    m03 = optimal(ind(3,3,"m0"));
    mp3 = optimal(ind(3,3,"mp"));
    ms3 = optimal(ind(3,3,"ms"));
    m03_unw = m03-param.mPL;
else
    % 1st stage
    m01 = optimal(ind(2,1,"m0"));
    mp1 = optimal(ind(2,1,"mp"));
    ms1 = optimal(ind(2,1,"ms"));
    m01_unw = m01-param.mPL;
    
    % 2nd stage
    m02 = optimal(ind(2,2,"m0"));
    mp2 = optimal(ind(2,2,"mp"));
    ms2 = optimal(ind(2,2,"ms"));
    m02_unw = m02-param.mPL;
end

if is_3stg
    % 3 stg input check - if total dv is delivered for std PL
    deliver_dv = g0*(Isp1*log(m01/(m01-mp1))+Isp2*log(m02/(m02-mp2))+Isp3*log(m03/(m03-mp3)));
else
    % 2 stg input check - if total dv is delivered for std PL
    deliver_dv = g0*(Isp1*log(m01/(m01-mp1))+Isp2*log(m02/(m02-mp2)));
end
if deliver_dv + 0.01 < param.vrq
    disp("ERROR - System with Insufficient dv"); return
end

if is_3stg
% iteratively solve for unmatch cases
dvf = @(mPL) g0*(Isp1*log((m01_unw+mPL)/(m01_unw+mPL-mp1))+ ...
    Isp2*log((m02_unw+mPL)/(m02_unw+mPL-mp2))+ ...
    Isp3*log((m03_unw+mPL)/(m03_unw+mPL-mp3)));
else
dvf = @(mPL) g0*(Isp1*log((m01_unw+mPL)/(m01_unw+mPL-mp1))+ ...
    Isp2*log((m02_unw+mPL)/(m02_unw+mPL-mp2)));
end

MPL_SOLVE = nan(size(OP_SITE_LAT));
for lat = 1:size(OP_SITE_LAT,2)
    for inc = size(OP_SITE_LAT,1):-1:lat
        if isnan(REQ_UNMATCH_DV_BUDGET(inc,lat))
            continue
        end
        ddv = @(mPL) dvf(mPL)-REQ_UNMATCH_DV_BUDGET(inc,lat);
        MPL_SOLVE(inc,lat) = fzero(ddv,param.mPL);
    end
end

if plot_imperial
    mpl = param.mPL*2.20462;
    MPL = MPL_SOLVE*2.20462;
    bound_pl_low = param.mPL_lowerBound*2.20462;
    bound_pl_high = param.mPL_highBound*2.20462;
    label_mpl = "Payload capacity (lb)";

    low_bound = param.mPL_lowerBound*2.20462;
    high_bound = param.mPL_highBound*2.20462;
else
    mpl = param.mPL;
    MPL = MPL_SOLVE;
    bound_pl_low = param.mPL_lowerBound;
    bound_pl_high = param.mPL_highBound;
    label_mpl = "Payload capacity (kg)";

    low_bound = param.mPL_lowerBound;
    high_bound = param.mPL_highBound;
end

figure; surf(OP_SITE_LAT,OP_ORB_INC,MPL,'EdgeColor','none'); hold on
surf(OP_SITE_LAT,OP_ORB_INC,ones(size(OP_SITE_LAT))*low_bound,'EdgeColor','none','FaceColor',[107 175 226]/255,'FaceAlpha',0.2);
surf(OP_SITE_LAT,OP_ORB_INC,ones(size(OP_SITE_LAT))*high_bound,'EdgeColor','none','FaceColor',[107 175 226]/255,'FaceAlpha',0.2);
xlabel("$\phi_{site}$ ($^{\circ}$)"); ylabel("$i_{orbit}$ ($^{\circ}$)");
zlabel(label_mpl);

figure; contourf(OP_SITE_LAT,OP_ORB_INC,MPL, 41, 'EdgeColor','none','FaceAlpha',1); hold on
% xline(30,"k--",'LineWidth',1.2); yline(100,"k--",'LineWidth',1.2); 
% scatter(30,100,'filled','MarkerFaceColor',[247 129 52]/255,'LineWidth',1.2);
clim([bound_pl_low bound_pl_high]); c = colorbar; c.Label.String = label_mpl;
xlabel("Launch site latitude ($^{\circ}$)"); ylabel("Desired orbit inclination ($^{\circ}$)");
title(title_name)

%% Helper function ind for 3 stages
function i = ind(total_stage,stage,mass_type)
    ind_m0 = 1; ind_ms = 2; ind_mp = 3; num_param = 5;
    if mass_type == "m0"
        i = (total_stage-stage)*num_param + ind_m0;
    elseif mass_type == "ms"
        i = (total_stage-stage)*num_param + ind_ms;
    elseif mass_type == "mp"
        i = (total_stage-stage)*num_param + ind_mp;
    else
        disp("ERROR - Not Recognizable Mass Type")
    end
end

end