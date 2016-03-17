% filename: Breakthough.m
% author: Marc Hesse
% date: 17 Mar 2016
% Description: Illustrate the dependence of the breakthough time on
%              input parameters.

clear all, close all, clc
%% Paramters we think we know
Param.mu_w = 8.90e-4; % Pas
Param.mu_g = 1.48e-5; % Pas
Param.s_wr = 0.35;    % [-]

%% Assume quadratic relative permeabilities
Param.n_g = 3;
Param.n_w = 1.25;
Nk = 60; Nr = 20; Ns = 15;

k_rg_end =linspace(0,1,Nk+1);         k_rg_end(1) = [];
s_gr = linspace(0,1-Param.s_wr,Nr+1); s_gr(end)     = []; 
[K,S] = meshgrid(k_rg_end,s_gr);
TB = K*0;

for j=1:Nk
    for i=1:Nr
        Param.k_rg_end = K(i,j);  
        Param.s_gr = S(i,j);
        [v,sg,vs,ss] = BuckleyLeverett(Param,Ns);
        TB(i,j) = 1/vs;
%         plotBL(Param)
%         pause(.1)
    end
end

figure
subplot 121
contour(K,S,TB), colorbar
xlabel('k_{rg}^*')
ylabel('s_{gr}')
xlim([0 1])
pbaspect([1 1 1])

%% Assume sgr = 0.2 and k_rg_end = 0.1
Param.k_rg_end = 0.05;  
Param.s_gr = 0.25;

Ng = 20; Nw = 20; 

n_g = linspace(1,4,Ng);       
n_w = linspace(1,4,Nw);
[G,W] = meshgrid(n_g,n_w);
TB = G*0;

for i=1:Ng
    for j=1:Nw
        Param.n_g = G(i,j);  
        Param.n_w = W(i,j);
        [v,sg,vs,ss] = BuckleyLeverett(Param,Ns);
        TB(i,j) = 1/vs;
%         plotBL(Param)
%         pause(.1)
    end
end

subplot 122
contour(G,W,TB), colorbar
xlabel('n_g')
ylabel('n_w')
pbaspect([1 1 1])
