Ns=10;
numE=2;
Egap=zeros(1,length(Ns));
for j=1:length(Ns)
    j
    J=ones(1,Ns(j));
    EB=zeros(1,Ns(j));
    H=hsb(J,EB,0);
    Erg=eigs(H,numE,'sa');
    Egap(j)=Erg(2)-Erg(1);
end
    

