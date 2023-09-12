clear; clc; close all

%% sizing iteration - two stage rocket
param = sysParam();
vrq_ratio_stg_1 = 0.01:0.001:0.99;
vrq_ratio_stg_2 = 1-vrq_ratio_stg_1;

dvrq_stg_1 = vrq_ratio_stg_1*param.vrq;
dvrq_stg_2 = vrq_ratio_stg_2*param.vrq;

sizing_table = nan([10 size(vrq_ratio_stg_1,2)]);
for i = 1:size(vrq_ratio_stg_1,2)
    sizing_table(:,i) = TWTOSizing(dvrq_stg_2(i),dvrq_stg_1(i),...
        param.mPL,param.sigma_stg1,param.sigma_stg2,param.Isp_stg1*param.g0,param.Isp_stg2*param.g0);
end

figure; subplot(2,1,1)
plot(vrq_ratio_stg_1,sizing_table(6,:));
subplot(2,1,2)
plot(vrq_ratio_stg_1,sizing_table(1,:)./sizing_table(6,:));

%% sizing iteration - three stage rocket
vrq_ratio_stg_1 = 0.01:0.005:0.99;
vrq_ratio_stg_2 = 0.01:0.005:0.99;
[VRQ1,VRQ2] = meshgrid(vrq_ratio_stg_1,vrq_ratio_stg_2);
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

% figure; contourf(VRQ1,VRQ2,sizing_table_3stg(:,:,11),100,LineStyle="none");
% colorbar; clim([2000,50000])
% title("Total Liftoff Weight"); xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
% figure; contourf(VRQ1,VRQ2,sizing_table_3stg(:,:,6),100,LineStyle="none");
% colorbar; clim([0,5000])
% title("Second Stage Weight"); colorbar; xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
% figure; contourf(VRQ1,VRQ2,sizing_table_3stg(:,:,1),100,LineStyle="none");
% colorbar; clim([0,1000])
% title("Third Stage Weight"); colorbar; xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");

figure; plot3(VRQ1,VRQ2,sizing_table_3stg(:,:,11)); zlim([0,1e4])
title("Total Liftoff Weight"); xlabel("% 1st stage dV"); ylabel("% 2nd stage dV");
