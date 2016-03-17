function []=plotBL(Param)
% name: plotBL.m
% author: Marc Hesse
% date: 17 Mar 2016

% %% Paramters we think we know
% Param.mu_w = 8.90e-4; % Pas
% Param.mu_g = 1.48e-5; % Pas
% Param.s_wr = 0.35;    % [-]
% 
% %% Assume quadratic relative permeabilities
% Param.n_g = 2;
% Param.n_w = 2;
% Ns = 10;
% 
% Param.k_rg_end =.1;   
% Param.s_gr = .2; 

[fg,fw] = FractionalFlow(Param);
[v,sg,vs,ss] = BuckleyLeverett(Param,1e2);

subplot(1,3,1)
s = linspace(0,1,1e2); 
plot(s,fg(s)), hold on
plot([0 ss],[0,fg(ss)],'k-'), hold off
xlim([0 1]), ylim([0 1])
pbaspect([1 .8 1])

subplot(1,3,2)
plot(v,sg)
ylim([0 1])
pbaspect([1 .8 1])

subplot(1,3,3)
plot(1./v,sg)
xlim([0 1]), ylim([0 1])
pbaspect([1 .8 1])
