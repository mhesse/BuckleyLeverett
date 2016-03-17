function [v,sg,vs,ss] = BuckleyLeverett(Param,N)
% name: BuckleyLeverett.m
% author: Marc Hesse
% date: 03/10/2016
% description:
% Plot BL profile to explore effect of coefficients
mu_w = Param.mu_w;
mu_g = Param.mu_g;
k_rg_end = Param.k_rg_end;
s_wr     = Param.s_wr;
s_gr     = Param.s_gr; 
n_g      = Param.n_g;
n_w      = Param.n_w;

k_rg = @(sg) k_rg_end*((sg-s_gr)./(1-s_wr-s_gr)).^n_g.*(sg<1-s_wr).*(sg>s_gr)+k_rg_end*(sg>1-s_wr)+0*(sg<s_gr);
k_rw = @(sw) ((sw-s_wr)/(1-s_wr)).^n_w.*(sw>s_wr)+0;

sg = linspace(0,1,1e3); sw = 1-sg;

lam_w = @(sw) k_rw(sw)./mu_w;
lam_g = @(sg) k_rg(sg)./mu_g;

fg = @(sg) lam_g(sg)./(lam_g(sg)+lam_w(1-sg));
fw = @(sw) lam_w(sw)./(lam_g(1-sw)+lam_w(sw));

% Rel. perm. derivatives
dkgdsg = @(sg) n_g*k_rg_end*(sg-s_gr).^(n_g-1)/(1-s_gr-s_wr).^n_g.*(sg<1-s_wr).*(sg>s_gr)+0*(sg>1-s_wr)+0*(sg<s_gr);
dkwdsw = @(sw) n_w*(sw-s_wr).^(n_w-1)/(1-s_wr).^n_w.*(sw>s_wr)+0;
dkwdsg = @(sg) -dkwdsw(1-sg);

% Mobility derivatives
dlamgdsg = @(sg) dkgdsg(sg)/mu_g;
dlamwdsw = @(sw) dkwdsw(sw)/mu_w;
dlamwdsg = @(sg) dkwdsg(sg)/mu_w;

% Fractional flow derivative
dfgdsg = @(sg) ( dlamgdsg(sg).*(lam_g(sg)+lam_w(1-sg)) - lam_g(sg).*(dlamgdsg(sg)+dlamwdsg(sg)) )./(lam_g(sg)+lam_w(1-sg)).^2;

%% Find tangent point
% Find max slope
sg_max = fminbnd(@(sg) -dfgdsg(sg),0,1);

% Find tangent point
tan_obj = @(sg,n) (dfgdsg(sg) - fg(sg)./sg).^n;  % Tangent point objective function
sg_tan_min = fminbnd(@(sg) tan_obj(sg,1),sg_max,1);
sg_tan = fminbnd(@(sg) tan_obj(sg,2),sg_max,sg_tan_min,optimset('TolX',1e-12,'Display','off'));
% fprintf('sg_s = %3.5f: Lambda = %3.5f dfds = %3.5f\n',sg_tan,fg(sg_tan)/sg_tan,dfgdsg(sg_tan))

%% Buckley-Leverett solution
sg_rar = linspace(sg_tan,1);
v_rar = dfgdsg(sg_rar);
vs = fg(sg_tan)/sg_tan;

%% Assemble output
sg = [0,sg_rar]';
v  = [vs,v_rar]';
ss = sg_tan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotting
% figure('name','Relative permeabilities')
% subplot(2,3,1)
% plot(sw,k_rw(sw),'r-'), hold on
% plot(s_wr,0,'ro','markerfacecolor','w')
% xlabel('s_{w}'), pbaspect([1 .8 1])
% ylabel('k_{rw}')
% 
% subplot(2,3,2)
% plot(sg,k_rg(sg),'b-'), hold on
% plot([s_gr,1-s_wr],[0 k_rg_end],'bo','markerfacecolor','w')
% xlabel('s_{g}'), ylim([0 1]), pbaspect([1 .8 1])
% ylabel('k_{rg}')
% 
% subplot(2,3,3)
% plot(sg,k_rg(sg),'b-'), hold on
% plot(sg,k_rw(sw),'r-')
% plot(1-s_wr,k_rg_end,'bo','markerfacecolor','w')
% plot(1-s_wr,0,'ro','markerfacecolor','w')
% xlabel('s_{g}'), ylim([0 1]), pbaspect([1 .8 1])
% legend('k_{rg}','k_{rw}')
% 
% subplot(2,3,4)
% sw_ave = (sw(2:end)+sw(1:end-1))/2; dsw = sw(2)-sw(1);
% plot(sw,dkwdsw(sw),'b-',sw_ave,diff(k_rw(sw))/dsw,'k--')
% xlabel('s_{w}'), pbaspect([1 .8 1])
% ylabel('dk_{rw}/ds_w')
% 
% subplot(2,3,5)
% sg_ave = (sg(2:end)+sg(1:end-1))/2; dsg = sg(2)-sg(1);
% plot(sg,dkgdsg(sg),'b-',sg_ave,diff(k_rg(sg))/dsg,'k--')
% xlabel('s_{g}'),  pbaspect([1 .8 1])
% ylabel('dk_{rg}/ds_g')
% 
% subplot(2,3,6)
% plot(sg,dkgdsg(sg),'b-'), hold on
% plot(sg,dkwdsg(sg),'r-')
% xlabel('s_{g}'), pbaspect([1 .8 1])
% legend('dk_{rg}/ds_g','dk_{rw}/ds_w')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('name','Fractional flow')
% subplot(2,3,1)
% plot(sg,sg,'k:',sg,fg(sg),'b-'), hold on
% plot([0 sg_tan],[0 fg(sg_tan)],'k-')
% xlabel 's_g', ylabel 'f_g'
% pbaspect([1 1 1])
% 
% subplot(2,3,2)
% plot(sw,sw,'k:',sw,fw(sw),'r-')
% xlabel 's_w', ylabel 'f_w'
% pbaspect([1 1 1])
% 
% subplot(2,3,3)
% plot(sg,tan_obj(sg,1),[0 1],[0 0],'k:'), hold on
% plot(sg_max,0,'o')
% plot(sg_tan_min,tan_obj(sg_tan_min,1),'o')
% plot(sg_tan,tan_obj(sg_tan,1),'k.','markersize',12)
% 
% subplot(2,3,4)
% plot(sg,dfgdsg(sg),'b-',sg_max,dfgdsg(sg_max),'bo'), hold on
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('name','Buckley-Leverett solution')
% subplot(1,3,1)
% plot(sg,sg,'k:',sg,fg(sg),'b-'), hold on
% plot([0 sg_tan],[0 fg(sg_tan)],'k-')
% plot(sg_tan,fg(sg_tan),'ko','markerfacecolor','w')
% xlabel 's_g', ylabel 'f_g'
% xlim([0 1]), ylim([0 1])
% pbaspect([1 1 1])
% 
% subplot(1,3,2)
% plot(v_rar,sg_rar,'k-'), hold on
% plot(v_s*[1 1],[0,sg_tan],'k-')
% xlim([0 3]), ylim([0 1])
% pbaspect([1 .8 1])
% xlabel 'v', ylabel 's_g'
% 
% subplot(1,3,3)
% plot(1./v_rar,sg_rar,'k-','markerfacecolor','w'), hold on
% plot(1/v_s*[1 1],[0,sg_tan],'k-')
% xlim([0 2]), ylim([0 1])
% set(gca,'xtick',[0:.2:2])
% pbaspect([1 .8 1])
% xlabel 'pv', ylabel 's_g'


