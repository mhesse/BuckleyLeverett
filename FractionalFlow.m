function [fg,fw] = FractionalFlow(Param)
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
fw = @(sg) lam_w(sw)./(lam_g(1-sw)+lam_w(sw));