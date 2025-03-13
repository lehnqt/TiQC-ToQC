function [iF,iG]=mps_infid_T_robust(J,H0,Hc,x,mps0,mpstg,sv_min,D,midstep,nt,iscpr,iso,varT)
T=x(end);
c=x(1:end-1);
[iF,iG1]=mps_infid_robust(J,H0,Hc,c,T,mps0,mpstg,sv_min,D,midstep,nt,iscpr,iso);
if varT==1
dt=10^(-10);    
T=T+dt;
iF2=mps_infid_nograd_robust(J,H0,Hc,c,T,mps0,mpstg,sv_min,D,midstep,nt,iscpr,iso,1);
iGT=(iF2-iF)/dt;
else
iGT=0;
end
iG=[iG1;iGT];
end
