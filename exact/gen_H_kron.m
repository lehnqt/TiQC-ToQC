function [H0,Hc]=gen_H_kron(H_cert,H_ctrl,sys_num,sys_dim)
cert_num=length(H_cert);
ctrl_num=length(H_ctrl);
H0=struct('op',cell(1,cert_num),'ft',cell(1,cert_num));
Hc=struct('op',cell(1,ctrl_num),'ft',cell(1,ctrl_num));
M=sys_dim;
for j=1:cert_num
H0(j).ft=H_cert(j).ft;
H0(j).op=sparse(prod(M),prod(M));
    for k=1:size(H_cert(j).op,1)
        H_temp=1;
        count=0;
        for l=1:sys_num
            if ismember(l,H_cert(j).sys)
                count=count+1;
            H_temp=kron(H_temp,H_cert(j).op{k,count});
            else
                H_temp=kron(H_temp,speye(M(l)));
            end
        end
        H0(j).op=H0(j).op+H_temp;
    end
        if ~issparse(H0(j).op)
        fprintf('warning: H0 %d is not sparse\n',j);
        end
end
for j=1:ctrl_num
Hc(j).ft=H_ctrl(j).ft;
Hc(j).op=sparse(prod(M),prod(M));
    for k=1:size(H_ctrl(j).op,1)
        H_temp=1;
        count=0;
        for l=1:sys_num
            if ismember(l,H_ctrl(j).sys)
                count=count+1;
            H_temp=kron(H_temp,H_ctrl(j).op{k,count});
            else
                H_temp=kron(H_temp,speye(M(l)));
            end
        end
        Hc(j).op=Hc(j).op+H_temp;
    end
    if ~issparse(Hc(j).op)
       fprintf('warning: Hc %d is not sparse\n',j);
    end
end
end