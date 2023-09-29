clear; clc; close all

%% edge case sizing - assume 1% structure SSTO with 24 kg PL
% propellant selection
Isp = 270;
density = 1370;

% payload mass and structure ratio
mPL = 24; sigma = 0.01; g0 = 9.80655;
% assume 92% of volume is propellant
Volup = (24*0.0254)^2*(240-15)*0.0254;
% propellant mass
mp = Volup*density;
% structure mass
ms = (mp+mPL)/(1-sigma)*sigma;

vfinal = g0*Isp*log((ms+mp+mPL)/(ms+mPL));
disp("dv available: "+vfinal)
