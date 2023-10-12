clear; clc; close all

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',16)

load("SatDB.mat")
% clean the data
incli = SatDB2023(:,5);
incli(incli>180) = nan;

figure; histogram(incli,"EdgeColor","auto","FaceColor",[0 104 56]/255,"FaceAlpha",1);
hold on; grid on; text(120, 3900, "Courtesy: Union of Concerned Scientists", "Color",[0.8 0.8 0.8])
xlabel("Orbit Inclination ($^{\circ}$)"); ylabel("Satellite Count");
title("On-Orbit Satellite Count on Orbit Inclindation")

sat_count = [sum(incli<101) sum(incli>=101)]; explode = [0,1];
figure; p = pie3(sat_count,explode);
p(3).FaceColor = [0 104 56]/255; 
p(2).FaceColor = [0 104 56]/255; p(2).EdgeColor = "none";
p(1).FaceColor = [0 104 56]/255; p(1).EdgeColor = "none";
p(6).FaceColor = [247 129 52]/255;
p(5).FaceColor = [247 129 52]/255; p(5).EdgeColor = "none";
p(7).FaceColor = [247 129 52]/255; p(7).EdgeColor = "none";
l = legend("6694 Satellites with inclination $\le$ $100^{\circ}$",...
    "21 Satellites with inclination $> 100^{\circ}$","interpreter","latex","Location","south");