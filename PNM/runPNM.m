%run PNM simulations
clc;clear; close all

outputPath = './Result/';
mkdir(outputPath)
load pn.mat
L = 0.0368; %packing size

%% rheology model PHPA06
rho = 1000.6;
mu_0 = 16.14;
mu_inf = 0.0078505;
m = 0.39127;
n = 0.39616;
cm = PowerLawModel(mu_0, mu_inf, m, n);

%% pore-network
pn = PoreNetwork(L, throatConn, throatProp, nodeProp, cm);
dP_max = rho*9.8*L*2;
dP_min = rho*9.8*L*1e-4;
num = 500;
maxIter = 100;
tol = 1e-6;
[Q, dP, ratio_stat2] = pn.SolvePQcurve(dP_min, dP_max, num, maxIter, tol, outputPath);

