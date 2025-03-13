Ns=20;
dJ=0;
numE=2;
EB=zeros(1,Ns);
Erg=zeros(length(dJ),numE);
for j=1:length(dJ)
    j
    J=-(1+dJ(j)*((-1).^(0:Ns-1)));
    H=hsb(J,EB,1);
    Erg(j,:)=eigs(H,numE,'sa');
end
    

