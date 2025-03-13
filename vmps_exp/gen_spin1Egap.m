Ns=8:2:14;
numE=6;
Erg=zeros(numE,length(Ns));
Egap=zeros(1,length(Ns));
for j=1:length(Ns)
    j
    J=ones(1,Ns(j));
    EB=zeros(1,Ns(j));
    H=spin1_hsb(J,EB,1);
    Erg(:,j)=eigs(H,numE,'sa');
    Egap(j)=Erg(2,j)-Erg(1,j);
end
    

