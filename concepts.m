clear; clc; close all

enable = [false, false, false, false, false, true];

%% Concept 1 - Thumper - solid solid
if enable(1)
disp("##### Concept 1 - Thumper - solid solid #####")
is_3stg_c1 = false;
param_c1 = sysParam();
param_c1.mPL = 25;
% first stage
param_c1.Isp_stg1 = 290; % CL-20 (near 100% effective Isp)
param_c1.sigma_stg1 = 0.06;
param_c1.density_stg1 = 2040;
% second stage
param_c1.Isp_stg2 = 290; % CL-20 (near 100% effective Isp)
param_c1.sigma_stg2 = 0.08;
param_c1.density_stg2 = 2040;
% ram/scram condition
param_c1.is_scram = false;
param_c1.is_scram_solid_boost = false;
% concept sizing
[optimal_c1, dvdisb_c1] = iterativeConceptDesign(param_c1,is_3stg_c1);
end

%% Concept 2 - Energizer - dual ram/scram RDRE
if enable(2)
disp("##### Concept 2 - Energizer - dual ram/scram RDRE #####")
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
end

%% Concept 3 - Bugz - solid RDRE
if enable(3)
disp("##### Concept 3 - Bugz - solid RDRE #####")
is_3stg_c3 = false;
param_c3 = sysParam();
param_c3.mPL = 40;
% first stage
param_c3.Isp_stg1 = 290; % CL-20 (near 100% effective Isp)
param_c3.sigma_stg1 = 0.06;
param_c3.density_stg1 = 2040;
% second stage
param_c3.Isp_stg2 = 365; % RP-1 H2O2 (90% effective Isp with RDRE)
param_c3.sigma_stg2 = 0.1;
param_c3.density_stg2 = 1370; % O/F ratio 7.07
% ram/scram condition
param_c3.is_scram = false;
param_c3.is_scram_solid_boost = false;
% concept sizing
[optimal_c3, dvdisb_c3] = iterativeConceptDesign(param_c3,is_3stg_c3);
end

%% Concept 4 - Easter - solid solid solid
if enable(4)
disp("##### Concept 4 - Easter - solid solid solid #####")
is_3stg_c4 = true;
param_c4 = sysParam();
param_c4.mPL = 35;
% first stage
param_c4.Isp_stg1 = 290; % CL-20 (near 100% effective Isp)
param_c4.sigma_stg1 = 0.06;
param_c4.density_stg1 = 2040;
% second stage
param_c4.Isp_stg2 = 290; % CL-20 (near 100% effective Isp)
param_c4.sigma_stg2 = 0.08;
param_c4.density_stg2 = 2040;
% third stage
param_c4.Isp_stg3 = 290; % CL-20 (near 100% effective Isp)
param_c4.sigma_stg3 = 0.1;
param_c4.density_stg3 = 2040;
% ram/scram condition
param_c4.is_scram = false;
param_c4.is_scram_solid_boost = false;
% concept sizing
[optimal_c4, dvdisb_c4] = iterativeConceptDesign(param_c4,is_3stg_c4);
end

%% Concept 5 - Peter - solid scramjet RDRE
if enable(5)
disp("##### Concept 5 - Peter - solid scramjet RDRE #####")
is_3stg_c5 = true;
param_c5 = sysParam();
param_c5.mPL = 45;
% first stage
param_c5.Isp_stg1 = 290; % CL-20 (near 100% effective Isp)
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
end

%% Concept 6 - XXX - solid solid RDRE
if enable(6)
disp("##### Concept 6 - XXX - solid solid RDRE #####")
is_3stg_c6 = true;
param_c6 = sysParam();
param_c6.mPL = 42;
% first stage
param_c6.Isp_stg1 = 265; % CL-20 (near 100% effective Isp)
param_c6.sigma_stg1 = 0.06;
param_c6.density_stg1 = 2040;
% second stage
param_c6.Isp_stg2 = 265; % CL-20 (near 100% effective Isp)
param_c6.sigma_stg2 = 0.06;
param_c6.density_stg2 = 2040;
% third stage
param_c6.Isp_stg3 = 365; % CL-20 (near 100% effective Isp)
param_c6.sigma_stg3 = 0.1;
param_c6.density_stg3 = 1370;
% ram/scram condition
param_c6.is_scram = false;
param_c6.is_scram_solid_boost = false;
% concept sizing
[optimal_c6, dvdisb_c6] = iterativeConceptDesign(param_c6,is_3stg_c6);
end
