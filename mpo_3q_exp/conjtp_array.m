function tns1=conjtp_array(tns)%calculate conjugate transpose of a tnoperator
[m,n]=size(tns);
%complex transpose of tns
for j=1:m
    for k=1:n
    tns{j,k}=permute(conj(tns{j,k}),[1,3,2,4]);
    end
end
tns1=tns;
end