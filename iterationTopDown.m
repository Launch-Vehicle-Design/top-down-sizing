clear; clc; close all

param = sysParam();
max_dv_stg_1_perc = 0.99;
if param.is_scram
    lambda = 1.401; R = 287;
    estimate_losses = param.vrq - param.orb_vel;
    max_dv_stg_1 = param.scram_mach*sqrt(lambda*R*param.scram_ceiling_temp) + 0.9*estimate_losses;
    max_dv_stg_1_perc = max_dv_stg_1/param.vrq;
end

solid_boost_dv_stg_1_perc = 0;
if param.is_scram && param.is_scram_solid_boost
    lambda = 1.4; R = 287;
    solid_boost_dv = param.scram_solid_boost_mach*sqrt(lambda*R*param.release_temp);
    solid_boost_dv_stg_1_perc = solid_boost_dv/param.vrq;
end

%% sizing iteration - two stage rocket
vrq_ratio_stg_1 = 0.01:0.001:max_dv_stg_1_perc;
vrq_ratio_stg_2 = 1-vrq_ratio_stg_1;

dvrq_stg_1 = vrq_ratio_stg_1*param.vrq;
dvrq_stg_2 = vrq_ratio_stg_2*param.vrq;

sizing_table = nan([10 size(vrq_ratio_stg_1,2)]);
for i = 1:size(vrq_ratio_stg_1,2)
    sizing_table(:,i) = TWTOSizing(dvrq_stg_2(i),dvrq_stg_1(i),...
        param.mPL,param.sigma_stg1,param.sigma_stg2,param.Isp_stg1*param.g0,param.Isp_stg2*param.g0);
end
fuel_vol_2stg = sizing_table(3,:)/param.density_stg2 + sizing_table(8,:)/param.density_stg1;

figure; subplot(2,1,1)
plot(vrq_ratio_stg_1,sizing_table(6,:));
xlabel("Stage 1 dv%"); ylabel("Total vehicle mass");
ylim([0 param.bounding_mass]); grid on
subplot(2,1,2)
plot(vrq_ratio_stg_1,sizing_table(1,:)./sizing_table(6,:));
xlabel("Stage 1 dv%"); ylabel("stage 2 mass / total mass"); grid on

%% sizing iteration - three stage rocket
vrq_ratio_stg_1_3stg = 0.01:0.002:0.99;
if param.is_scram_solid_boost
    vrq_ratio_stg_1_3stg = solid_boost_dv_stg_1_perc;
end
vrq_ratio_stg_2_3stg = 0.01:0.002:max_dv_stg_1_perc-solid_boost_dv_stg_1_perc;
[VRQ1,VRQ2] = meshgrid(vrq_ratio_stg_1_3stg,vrq_ratio_stg_2_3stg);
VRQ3 = 1-VRQ1-VRQ2; VRQ3(VRQ3<0) = nan;

dVRQ1 = VRQ1*param.vrq; dVRQ2 = VRQ2*param.vrq;
dVRQ3 = VRQ3*param.vrq;

sizing_table_3stg = nan([size(dVRQ1,1) size(dVRQ1,2) 15]);
for i = 1:size(dVRQ1,1)
    for j = 1:size(dVRQ1,2)
        sizing_vec = TWTOSizing(dVRQ3(i,j),dVRQ2(i,j), ...
            param.mPL,param.sigma_stg2,param.sigma_stg3,param.Isp_stg2*param.g0,param.Isp_stg3*param.g0);
        % first stage
        mPL_stg_1 = sizing_vec(6); sigma_stg_1 = param.sigma_stg1; ve_stg_1 = param.Isp_stg1*param.g0;
        [m0_stg_1,ms_stg_1,mp_stg_1,mu_stg_1,PI_stg_1] = sizing(dVRQ1(i,j),mPL_stg_1,sigma_stg_1,ve_stg_1);
        sizing_vec(11:15) = [m0_stg_1 ms_stg_1 mp_stg_1 mu_stg_1 PI_stg_1]';
        if sum(sizing_vec < 0) ~= 0
            sizing_vec = nan(15,1);
        end
        sizing_table_3stg(i,j,:) = sizing_vec;
    end
end

if ~(param.is_scram && param.is_scram_solid_boost)
    figure; surf(VRQ1*100,VRQ2*100,sizing_table_3stg(:,:,11),"EdgeColor","interp");
    clim([0,param.bounding_mass]); zlim([0,param.bounding_mass]); xlim([0,100]); ylim([0,100])
    title("Total Liftoff Weight"); xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
    
    fuel_vol_3stg = sizing_table_3stg(:,:,3)/param.density_stg3 + ...
        sizing_table_3stg(:,:,8)/param.density_stg2 + sizing_table_3stg(:,:,13)/param.density_stg1;
    figure; surf(VRQ1*100,VRQ2*100,fuel_vol_3stg,"EdgeColor","interp"); 
    clim([0,param.bounding_box_volu]); zlim([0,param.bounding_box_volu]); xlim([0,100]); ylim([0,100])
    title("Total Fuel Volume"); xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
    
    figure; surf(VRQ1*100,VRQ2*100,100*sizing_table_3stg(:,:,11)/param.bounding_mass,"EdgeColor","interp"); hold on
    surf(VRQ1*100,VRQ2*100,100*fuel_vol_3stg/param.bounding_box_volu,"EdgeColor","interp"); 
    clim([0,100]); zlim([0,100]); xlim([0,100]); ylim([0,100])
    xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
else
    figure; plot3(VRQ1*100,VRQ2*100,sizing_table_3stg(:,:,11));
    clim([0,param.bounding_mass]); zlim([0,param.bounding_mass]); xlim([0,100]); ylim([0,100])
    title("Total Liftoff Weight"); xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
    
    fuel_vol_3stg = sizing_table_3stg(:,:,3)/param.density_stg3 + ...
        sizing_table_3stg(:,:,8)/param.density_stg2 + sizing_table_3stg(:,:,13)/param.density_stg1;
    figure; plot3(VRQ1*100,VRQ2*100,fuel_vol_3stg); 
    clim([0,param.bounding_box_volu]); zlim([0,param.bounding_box_volu]); xlim([0,100]); ylim([0,100])
    title("Total Fuel Volume"); xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
end

close all
% locate minimum
if ~(param.is_scram && param.is_scram_solid_boost)
    [min_mass_2stg,ind_2stg] = min(sizing_table(6,:));
    fuel_vol_min_mass_2stg = fuel_vol_2stg(ind_2stg);
    disp("2 Stage - Min Mass = "+min_mass_2stg+" kg Fuel Vol = "+fuel_vol_min_mass_2stg+" m^3 as "+ ...
        fuel_vol_min_mass_2stg/param.bounding_box_volu*100+"% of the bounding box")
end

[min_mass_3stg,ind_3stg] = min(sizing_table_3stg(:,:,11),[],"all");
fuel_vol_min_mass_3stg = fuel_vol_3stg(ind_3stg);
disp("3 Stage - Min Mass = "+min_mass_3stg+" kg Fuel Vol = "+fuel_vol_min_mass_3stg+" m^3 as "+ ...
    fuel_vol_min_mass_3stg/param.bounding_box_volu*100+"% of the bounding box")
optimal_3stg = nan([size(sizing_table_3stg,3),1]);
for i = 1:size(sizing_table_3stg,3)
    temp = sizing_table_3stg(:,:,i);
    optimal_3stg(i) = temp(ind_3stg);
end
save("OptimalSolution.mat","optimal_3stg","param")