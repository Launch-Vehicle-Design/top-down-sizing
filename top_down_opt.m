clear; clc; close all
% optimization top down method
param.sigma = [0.08; 0.096];
param.c = [278; 377.2]*9.80665;
param.total_dv = 30000*0.3048;
param.mpl = 46;

cost = @(vs) log_cost_func(vs,param);
Aeq = ones(1,2); Beq = param.total_dv;
init_vs = [0.5; 0.5]*param.total_dv;

x = fmincon(cost,init_vs,[],[],Aeq,Beq,[],[]);

function val = cost_func(vs,param)
    val = param.mpl;
    for i = 1:length(vs)
        val = val/(exp(-vs(i)/param.c(i))-param.sigma(i))*(1-param.sigma(i));
    end
end

function val = log_cost_func(vs,param)
    val = 0;
    for i = 1:length(vs)
        val = val - log(exp(-vs(i)/param.c(i))-param.sigma(i));
    end
end