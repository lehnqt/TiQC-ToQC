function [iF,iG]=infid_rcp_T(H0,Hc,x,mpo0,mpotg,sv_min,D,midstep,nt,varT)
T=x(end);
c=x(1:end-1);
[iF,iG1]=infid_rcp(H0,Hc,c,T,mpo0,mpotg,sv_min,D,midstep,nt);
if varT==1
dt=10^(-10);    
T=T+dt;
iF2=infid_nograd_rcp(H0,Hc,c,T,mpo0,mpotg,sv_min,D,midstep,nt);
iGT=(iF2-iF)/dt;
else
iGT=0;
end
iG=[iG1;iGT];
end
