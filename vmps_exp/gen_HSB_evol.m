N=20
D=5
precision=1e-5;
dt=0.005;
A=2;
% jflipped=5;
% magnetization in z-direction

sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2);
s_plus=(sx+1i*sy)/2;
s_minus=(sx-1i*sy)/2;
p0=(id+sz)/2;
p1=(id-sz)/2;
oset1=cell(1,N);
for j=1:N
    oset1{1,j}=p0; 
end
oset2=cell(1,N);
for j=1:N
    oset2{1,j}=p1; 
end
oset3=cell(1,N);
for j=1:N
    oset3{1,j}=s_plus; 
end
oset4=cell(1,N);
for j=1:N
    oset4{1,j}=s_minus; 
end
% oset{1,jflipped}=sx;
% time evolution operator
h=kron(s_plus,s_minus)+kron(s_minus,s_plus)+A*(kron(id,sx))/2+A*(kron(sx,id))/2;
w=expm(-1i*dt*h);
w=reshape(w,[2,2,2,2]); w=permute(w,[1,3,2,4]); w=reshape(w,[4,4]);
[U,S,V]=svd2(w); eta=size(S,1);
U=U*sqrt(S); V=sqrt(S)*V;
U=reshape(U,[2,2,eta]); U=permute(U,[4,3,2,1]);
V=reshape(V,[eta,2,2]); V=permute(V,[1,4,3,2]);
I=reshape(id,[1,1,2,2]);
mpo_even=cell(1,N);
mpo_odd=cell(1,N);
for j=1:N
    mpo_even{j}=I; 
    mpo_odd{j}=I;
end
for j=1:2:(N-1) 
    mpo_odd{j}=U;
    mpo_odd{j+1}=V; 
end
for j=2:2:(N-1)
    mpo_even{j}=U;
    mpo_even{j+1}=V; 
end
% starting state (one spin flipped)
mps0=cell(1,N);
for j=1:N
% if j==jflipped 
%     state=[0; 1]; 
% else
    state=[1; 0];
% end
mps0{j}=reshape(state,[1,1,2]);
end
% time evolution
mps=mps0;
mzvalues=[];
Nstep=50;
for step=1:Nstep
fprintf('Step %2d:',step);
mps=reduceD(mps,mpo_even,D,precision);
mps=reduceD(mps,mpo_odd,D,precision);
mz=expectationvalue(mps,oset1)+expectationvalue(mps,oset2)+expectationvalue(mps,oset3)+expectationvalue(mps,oset4);
mzvalues=[mzvalues,mz];
fprintf('mz=%g\n',mz);
end
plot(1:Nstep,real(mzvalues));
