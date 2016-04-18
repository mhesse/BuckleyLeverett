% file name: plot_experiments.m
% author: Marc Hesse
% date: 17 Mar 2016
close all, clear all, clc
%load column_N2_01_25_2016
load column_N2_01_26_2016

% Pore volume
L = 123.7;   % [cm] column length 
R = 1/2; % [cm] column radius
A = pi*R^2; % [cm^2] cross-sectional area
V = L*A;     % [cm^3/ml] column volume
Vp = 40.8;   % [cm^3/ml] pore volume (email Toti: Wednesday, February 24, 2016 at 5:27 PM)
Q = 0.25;    % [cm^3/min] Volumetric injection rate
rhow = 1;    % [g/cm^3] Density of water

%% Derived quantities
phi = Vp/V;  % [-] Porosity of the column
tc = Vp/Q;   % [min] time to inject one pore volume
Mw = Vp*rhow;
time_ave = (time(1:end-1)+time(2:end))/2;
s_wr = 1-max(mass/Mw);

dmdt = diff(mass)./diff(time);

% Determine breakthrough
[dmdt_max,i_max]= max(dmdt);
tb = time_ave(i_max+1);
mass_tb = (mass(i_max+1)+mass(i_max+2))/2;

fprintf('phi = %3.2f s_wr = %3.2f\n',phi,s_wr)
fprintf('tb = %3.2f min = %3.2f PVI\n',tb,tb/tc)

%% Analytic solution
Param.mu_w = 8.90e-4; % Pas
Param.mu_g = 1.48e-5; % Pas
Param.s_wr = 0.35;    % [-]
Param.n_g = 3;
Param.n_w = 1.25;
Param.k_rg_end = .1;  
Param.s_gr = .25;
[mass_model,t_model]=CummulativeWater(Param,1e2);

subplot 211
plot(time/tc,mass/Mw), hold on
plot(t_model,mass_model)
plot(tb/tc,mass_tb/Mw,'k.')
plot([0 1],[0 1],'k:')
xlabel 'pore volume injected'
ylabel 'cum. mass/total mass'
xlim([0 8])
legend('data','label')
subplot 211


% plot(time,pres)%, ylim([0 .4])
xlabel 'time'
ylabel 'Pressure'
xlim([0 8])

figure
plotBL(Param)