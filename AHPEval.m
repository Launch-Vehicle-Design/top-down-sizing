function [PL_max_g, ms2mf, llf] = AHPEval(optimal,is_3stg,MPL,inc)
% assumed an initial thrust to weight ratio
T2W0 = [1.6 0.8 0.8]; g = 9.80665;
% max payload g loading
PL_max_g = 0; total_ms = 0; ms2mf = 1;
total_mpl_flex = 0;

total_stage = 2; 
if is_3stg
    total_stage = 3;
end
for i = 1:total_stage
    m0i = optimal(ind(total_stage,i,"m0"));
    mpi = optimal(ind(total_stage,i,"mp"));
    mfi = m0i-mpi; T_i = m0i*g*T2W0(i);
    g_si = T_i/mfi/g;
    if PL_max_g < g_si
        PL_max_g = g_si;
    end
    msi = optimal(ind(total_stage,i,"ms"));
    if i == total_stage
        ms2mf = total_ms/(total_ms+mfi);
    else
        total_ms = total_ms + msi;
    end
end

mpl = MPL; % mpl(isnan(mpl)) = 0;
load("SatDBTJ.mat", "Inclination", "SemiMajorAxis");
incli = Inclination(Inclination<180 & SemiMajorAxis<6378100+500000);
[sat_count,~] = histcounts(incli,0:180);
[data_count,~] = histcounts(inc,0:180);
ratio_sat_count = sat_count/length(incli);
start_ind = 1; llf = 0;
for i = 1:size(data_count,2)
    inds = start_ind:start_ind+data_count(i)-1;
    start_ind = start_ind + data_count(i);
    mpl_subset = mpl(inds,:); nan_purges_mpl = mpl_subset;
    nan_purges_mpl(isnan(nan_purges_mpl)) = 0;
    nan_purges_mpl(nan_purges_mpl > 100) = 100;
    llf = llf + sum(nan_purges_mpl,"all")/sum(~isnan(mpl_subset),"all")*ratio_sat_count(i);
end

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