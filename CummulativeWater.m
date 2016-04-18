function [mass,t]=CummulativeWater(Param,N)

[v,sg,vs,ss] = BuckleyLeverett(Param,N);
t = 1./v; 
t = [0;t];
sg = [0;sg];
mass = cumtrapz(t,1-sg);

tb = 1/vs;
