function [iF,iG]=infid_T(H0,Hc,x,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso,varT)
T=x(end);
c=x(1:end-1);
if nargout<2
   iF=infid_nograd(H0,Hc,c,T,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
else
[iF,iG1]=infid(H0,Hc,c,T,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
if varT==1
dt=10^(-10);    
T=T+dt;
iF2=infid_nograd(H0,Hc,c,T,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
iGT=(iF2-iF)/dt;
else
iGT=0;
end
iG=[iG1;iGT];
end
end
