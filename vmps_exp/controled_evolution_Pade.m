clear all;
N=7;
M=2;
A=1;
T=pi/2;
nbin=50;
dt=T/nbin;
J=ones(1,N-1);
% J(N-1)=0.01;
time_grid=dt*(1:nbin);
c=zeros(nbin,N,2);
for j=1:N
    c(:,j,1)=j*exp(-(time_grid-T/2).^2/(T/6));
%     ctemp=diff(c(:,j,1));
%     c(:,j,2)=[ctemp(1);ctemp] ;
    c(:,j,2)=c(:,j,1);
end
numK=25;
% c=ones(nbin,N,2);
H1kron =@(H1,j) kron(speye(M^(j-1)),kron(H1,speye(M^(N-j))));
H2kron =@(H1,H2,j,k) kron(speye(M^(j-1)),kron(kron(H1,speye(M^(k-j-1))),kron(H2,speye(M^(N-k)))));
% jflipped=5;
% magnetization in z-directionc=ones(nbin,N,2);
sx=sparse([0,1;1,0]); sy=sparse([0,-1i;1i,0]); sz=sparse([1,0;0,-1]); id=speye(2);
s_plus=(sx+1i*sy)/2;
s_minus=(sx-1i*sy)/2;
p0=(id+sz)/2;
p1=(id-sz)/2;
oset1=cell(1,N);
factor=@(j) (eq(j,1)+eq(j,N)+1)/2; 
O1=p0;
for j=1:N-1
    O1=kron(O1,p0); 
end
O2=p1;
for j=1:N-1
    O2=kron(O2,p1); 
end
O3=s_plus;
for j=1:N-1
    O3=kron(O3,s_plus); 
end
O4=s_minus;
for j=1:N-1
    O4=kron(O4,s_minus); 
end
OGHZ=(O1+O2+O3+O4)/2;
% oset{1,jflipped}=sx;
% time evolution operator
q0=[1;0];
q1=[0;1];
psi0=q0;
psi1=q1;
for j=1:N-1
    psi0=kron(psi0,q0);
    psi1=kron(psi1,q1);
end
% time evolution
psiGHZ=(psi0+psi1)/sqrt(2);
psi=psi0;
ghzoverlap_Pade=zeros(nbin,1);
H0=sparse(M^N,M^N);
for j=1:(N-1) 
%     h=kron(s_plus,s_minus)+kron(s_minus,s_plus)+(c(k,j,1)*(kron(sx,id))+c(k,j,2)*(kron(sy,id)))*factor(j)+(c(k,j+1,1)*kron(id,sx)+c(k,j+1,2)*kron(id,sy))*factor(j+1);
    H0=H0+J(j)*(H2kron(s_plus,s_minus,j,j+1)+H2kron(s_minus,s_plus,j,j+1));
end
Hc=cell(N,2);
for j=1:N-1 
%     h=kron(s_plus,s_minus)+kron(s_minus,s_plus)+(c(k,j,1)*(kron(sx,id))+c(k,j,2)*(kron(sy,id)))*factor(j)+(c(k,j+1,1)*kron(id,sx)+c(k,j+1,2)*kron(id,sy))*factor(j+1);
    Hc{j,1}=H1kron(sx,j);
    Hc{j,2}=H1kron(sy,j);
end
for k=1:nbin
    fprintf('step %d\n',k);
    Htot=H0;
for j=1:N-1
%     h=kron(s_plus,s_minus)+kron(s_minus,s_plus)+(c(k,j,1)*(kron(sx,id))+c(k,j,2)*(kron(sy,id)))*factor(j)+(c(k,j+1,1)*kron(id,sx)+c(k,j+1,2)*kron(id,sy))*factor(j+1);
    Htot=Htot+A*(c(k,j,1)*Hc{j,1}+c(k,j,2)*Hc{j,2});
end
psi=expm(-1i*dt*Htot)*psi;
% ghzoverlap_Pade(k)=psi'*OGHZ*psi;
ghzoverlap_Pade(k)=psiGHZ'*psi;
end
% starting state (one spin flipped)
plot(time_grid,c(:,1,1));
hold on;
plot(time_grid,c(:,1,2));
hold off;
figure;
plot(time_grid,real(ghzoverlap_Pade));
hold on;
plot(time_grid,imag(ghzoverlap_Pade));
hold on;
plot(time_grid,abs(ghzoverlap_Pade).^2);