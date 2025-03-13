function [iF,iG] = infid_T_fd(H0,Hc,x,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso,varT)
iF=infid_T_nograd(H0,Hc,x,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
if nargout>1
dx=10^(-10);
iG=zeros(length(x),1);
if varT==1
parfor jx=1:length(x)
    x1=x;
    x1(jx)=x1(jx)+dx;
    iF1=infid_T_nograd(H0,Hc,x1,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
    iG(jx)=(iF1-iF)/dx;
end
else
parfor jx=1:length(x)-1
    x1=x;
    x1(jx)=x1(jx)+dx;
    iF1=infid_T_nograd(H0,Hc,x1,mpo0,mpotg,sv_min,D,Dc,nsweep,midstep,nt,iscpr,iso);
    iG(jx)=(iF1-iF)/dx;
end
end
end
end