N_list=2:18;
numE=2;
E0_exact=zeros(1,length(N_list));
E1_exact=zeros(1,length(N_list));
fprintf('exact diagonalisation\n');
for j=1:length(N_list)
    N=N_list(j);
    fprintf('N=%d\n',N);
    J=ones(1,N);
    EB=zeros(1,N);
    H=hsb(J,EB,0);
    Erg=eigs(H,numE,'sa');
    E0_exact(j)=Erg(1);
    E1_exact(j)=Erg(2);
end
N_list2=2:30;
E0_mps=zeros(1,length(N_list2));
E1_mps=zeros(1,length(N_list2));
fprintf('mps calculation\n');
D=5;
precision=1e-5;
for k=1:length(N_list2)
   N=N_list2(k);
   fprintf('N=%d\n',N);
% Heisenberg Hamiltonian
M=3*(N-1);
hset=cell(M,N);
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2);
for m=1:M, for j=1:N, hset{m,j}=id; end; end
for j=1:(N-1)
hset{3*(j-1)+1,j}=sx; hset{3*(j-1)+1,j+1}=sx;
hset{3*(j-1)+2,j}=sy; hset{3*(j-1)+2,j+1}=sy;
hset{3*(j-1)+3,j}=sz; hset{3*(j-1)+3,j+1}=sz;
end
% ground state energy
randn('state',0);
[E0,mps0]=minimizeE(hset,D,precision,[]);
E0_mps(k)=E0;
% first excited state
[E1,mps1]=minimizeE(hset,D,precision,mps0);
E1_mps(k)=E1;
end