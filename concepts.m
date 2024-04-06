clear; clc; close all

enable = [false, false, true, false, false, false, false, false];

%% Concept 1 - Hopps - solid solid
if enable(1)
disp("##### Concept 1 - Hopps - solid solid #####")
is_3stg_c1 = false;
param_c1 = sysParam();
param_c1.mPL = 46;
% first stage
param_c1.Isp_stg1 = 280;
param_c1.sigma_stg1 = 0.085;
param_c1.density_stg1 = 1783.8;
% second stage
param_c1.Isp_stg2 = 320;
param_c1.sigma_stg2 = 0.096;
param_c1.density_stg2 = 1440;
% ram/scram condition
param_c1.is_scram = false;
param_c1.is_scram_solid_boost = false;
% concept sizing
[optimal_c1, dvdisb_c1] = iterativeConceptDesign(param_c1,is_3stg_c1);
[mpl_c1,~,incli_c1] = launchAziPLCapacity(optimal_c1, param_c1, is_3stg_c1, "Hopps");
[PL_max_g_c1, ms2mf_c1, llf_c1] = AHPEval(optimal_c1,is_3stg_c1,mpl_c1,incli_c1)
end

%% Concept 2 - Bugs - dual ram/scram RDRE
if enable(2)
disp("##### Concept 2 - Bugs - dual ram/scram RDRE #####")
is_3stg_c2 = false;
param_c2 = sysParam();
param_c2.mPL = 45;
% first stage
param_c2.Isp_stg1 = 1000; % JP-7 Ram/Scram avg effective Isp
param_c2.sigma_stg1 = 0.5; % Ram/scram structure ratio + no oxidier
param_c2.density_stg1 = 800; % JP-7
% second stage
param_c2.Isp_stg2 = 365; % RP-1 H2O2 (90% effective Isp with RDRE)
param_c2.sigma_stg2 = 0.1;
param_c2.density_stg2 = 1370; % O/F ratio 7.07
% ram/scram condition
param_c2.is_scram = true;
param_c2.is_scram_solid_boost = false;
param_c2.scram_ceiling = 75000;
param_c2.scram_ceiling_temp = 206;
param_c2.scram_mach = 9;
% concept sizing
[optimal_c2, dvdisb_c2] = iterativeConceptDesign(param_c2,is_3stg_c2);
[mpl_c2,~,incli_c2] = launchAziPLCapacity(optimal_c2, param_c2, is_3stg_c2, "Bugs");
[PL_max_g_c2, ms2mf_c2, llf_c2] = AHPEval(optimal_c2,is_3stg_c2,mpl_c2,incli_c2)
end

%% Concept 3 - Thumper - solid RDRE
if enable(3)
disp("##### Concept 3 - Thumper - solid RDRE #####")
is_3stg_c3 = false;
param_c3 = sysParam();
param_c3.mPL = 46;
% first stage
param_c3.Isp_stg1 = 278;
param_c3.sigma_stg1 = 0.085;
param_c3.density_stg1 = 1783.8;
% second stage
param_c3.Isp_stg2 = 369.5;
param_c3.sigma_stg2 = 0.096;
param_c3.density_stg2 = 1280;
% ram/scram condition
param_c3.is_scram = false;
param_c3.is_scram_solid_boost = false;
% concept sizing
[optimal_c3, dvdisb_c3] = iterativeConceptDesign(param_c3,is_3stg_c3);
[mpl_c3,~,incli_c3] = launchAziPLCapacity(optimal_c3, param_c3, is_3stg_c3, "Thumper");
[PL_max_g_c3, ms2mf_c3, llf_c3] = AHPEval(optimal_c3,is_3stg_c3,mpl_c3,incli_c3)

optimal = optimal_c3; dvdisb = dvdisb_c3; param = param_c3;
save("../lv-dynamics-model/trajectory_optimization/thumper_trajopt_selfcontain/thumper.mat",...
    "optimal","dvdisb","param");
end

%% Concept 4 - Trix - solid solid solid
if enable(4)
disp("##### Concept 4 - Trix - solid solid solid #####")
is_3stg_c4 = true;
param_c4 = sysParam();
param_c4.mPL = 22;
% first stage
param_c4.Isp_stg1 = 267.6; % CL-20 (near 100% effective Isp)
param_c4.sigma_stg1 = 0.06;
param_c4.density_stg1 = 2040;
% second stage
param_c4.Isp_stg2 = 267.6; % CL-20 (near 100% effective Isp)
param_c4.sigma_stg2 = 0.08;
param_c4.density_stg2 = 2040;
% third stage
param_c4.Isp_stg3 = 267.6; % CL-20 (near 100% effective Isp)
param_c4.sigma_stg3 = 0.1;
param_c4.density_stg3 = 2040;
% ram/scram condition
param_c4.is_scram = false;
param_c4.is_scram_solid_boost = false;
% concept sizing
[optimal_c4, dvdisb_c4] = iterativeConceptDesign(param_c4,is_3stg_c4);
[mpl_c4,~,incli_c4] = launchAziPLCapacity(optimal_c4, param_c4, is_3stg_c4, "Trix");
[PL_max_g_c4, ms2mf_c4, llf_c4] = AHPEval(optimal_c4,is_3stg_c4,mpl_c4,incli_c4)
end

%% Concept 5 - Energizer - solid scramjet RDRE
if enable(5)
disp("##### Concept 5 - Energizer - solid scramjet RDRE #####")
is_3stg_c5 = true;
param_c5 = sysParam();
param_c5.mPL = 45;
% first stage
param_c5.Isp_stg1 = 267.6; % CL-20 (near 100% effective Isp)
param_c5.sigma_stg1 = 0.06;
param_c5.density_stg1 = 2040;
% second stage
param_c5.Isp_stg2 = 1000; % CL-20 (near 100% effective Isp)
param_c5.sigma_stg2 = 0.5;
param_c5.density_stg2 = 800;
% third stage
param_c5.Isp_stg3 = 365; % CL-20 (near 100% effective Isp)
param_c5.sigma_stg3 = 0.1;
param_c5.density_stg3 = 1370;
% ram/scram condition
param_c5.is_scram = true;
param_c5.is_scram_solid_boost = true;
param.scram_solid_boost_mach = 5;
param_c5.scram_ceiling = 75000;
param_c5.scram_ceiling_temp = 206;
param_c5.scram_mach = 9;
% concept sizing
[optimal_c5, dvdisb_c5] = iterativeConceptDesign(param_c5,is_3stg_c5);
[mpl_c5,~,incli_c5] = launchAziPLCapacity(optimal_c5, param_c5, is_3stg_c5, "Energizer");
[PL_max_g_c5, ms2mf_c5, llf_c5] = AHPEval(optimal_c5,is_3stg_c5,mpl_c5,incli_c5)
end

%% Concept 6 - X27 - solid solid RDRE
if enable(6)
disp("##### Concept 6 - X27 - solid solid RDRE #####")
is_3stg_c6 = true;
param_c6 = sysParam();
param_c6.mPL = 45;
% first stage
param_c6.Isp_stg1 = 276; % CL-20 (near 95% effective Isp)
param_c6.sigma_stg1 = 0.06;
param_c6.density_stg1 = 2040;
% second stage
param_c6.Isp_stg2 = 276; % CL-20 (near 95% effective Isp)
param_c6.sigma_stg2 = 0.08;
param_c6.density_stg2 = 2040;
% third stage
param_c6.Isp_stg3 = 369.5;
param_c6.sigma_stg3 = 0.1;
param_c6.density_stg3 = 1370;
% ram/scram condition
param_c6.is_scram = false;
param_c6.is_scram_solid_boost = false;
% concept sizing
[optimal_c6, dvdisb_c6] = iterativeConceptDesign(param_c6,is_3stg_c6);
[mpl_c6,~,incli_c6] = launchAziPLCapacity(optimal_c6, param_c6, is_3stg_c6, "X27");
[PL_max_g_c6, ms2mf_c6, llf_c6] = AHPEval(optimal_c6,is_3stg_c6,mpl_c6,incli_c6)
end

%% Concept 7 - X28 - solid RDRE RDRE
if enable(7)
disp("##### Concept 7 - X28 - solid RDRE RDRE #####")
is_3stg_c7 = true;
param_c7 = sysParam();
param_c7.mPL = 45;
% first stage
param_c7.Isp_stg1 = 267.6; % CL-20 (near 95% effective Isp)
param_c7.sigma_stg1 = 0.06;
param_c7.density_stg1 = 2040;
% second stage
param_c7.Isp_stg2 = 369.5;
param_c7.sigma_stg2 = 0.1;
param_c7.density_stg2 = 1370;
% third stage
param_c7.Isp_stg3 = 369.5;
param_c7.sigma_stg3 = 0.12;
param_c7.density_stg3 = 1370;
% ram/scram condition
param_c7.is_scram = false;
param_c7.is_scram_solid_boost = false;
% concept sizing
[optimal_c7, dvdisb_c7] = iterativeConceptDesign(param_c7,is_3stg_c7);
[mpl_c7,~,incli_c7] = launchAziPLCapacity(optimal_c7, param_c7, is_3stg_c7, "X28");
[PL_max_g_c7, ms2mf_c7, llf_c7] = AHPEval(optimal_c7,is_3stg_c7,mpl_c7,incli_c7)
end

%% Concept 8 - X29 - RDRE RDRE
if enable(8)
disp("##### Concept 8 - X29 - RDRE RDRE #####")
is_3stg_c8 = false;
param_c8 = sysParam();
param_c8.mPL = 45;
% first stage
param_c8.Isp_stg1 = 369.5;
param_c8.sigma_stg1 = 0.1;
param_c8.density_stg1 = 1370;
% second stage
param_c8.Isp_stg2 = 369.5;
param_c8.sigma_stg2 = 0.12;
param_c8.density_stg2 = 1370;
% ram/scram condition
param_c8.is_scram = false;
param_c8.is_scram_solid_boost = false;
% concept sizing
[optimal_c8, dvdisb_c8] = iterativeConceptDesign(param_c8,is_3stg_c8);
[mpl_c8,~,incli_c8] = launchAziPLCapacity(optimal_c8, param_c8, is_3stg_c8, "X29");
[PL_max_g_c8, ms2mf_c8, llf_c8] = AHPEval(optimal_c8,is_3stg_c8,mpl_c8,incli_c8)
end
