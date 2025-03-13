t=(0:nbin-1)/nbin;
col=linspecer(N);
for k=1:N
plot(t,x_optmM(k,:)/Amax,'-o','MarkerFaceColor',col(k,:),'MarkerEdgeColor',col(k,:),'MarkerSize',4);
hold on;
end
hold off;
set(gca,'fontsize',14);
box on;
set(gcf,'color','w');
xlim([-0.05, 1]);
ylim([-1.05,1.05]);
legend(sprintfc('Q%d',1:N));

figure;
for k=1:N
plot(t,x_optmM(N+k,:)/Amax,'-o','MarkerFaceColor',col(k,:),'MarkerEdgeColor',col(k,:),'MarkerSize',4);
hold on;
end
hold off;
set(gca,'fontsize',14);
box on;
set(gcf,'color','w');
xlim([-0.05,1]);
ylim([-1.05,1.05]);
legend(sprintfc('Q%d',1:N));

F_optm=1-fun(x_optm);
[~,i_optm]=min(F_optm);
P=zeros(M,nbin+1,N);
for k=1:N
    for j=1:M
        P(j,1:nbin+1,k)=pop_evol(H0,Hc,Delta_t,proj,Pq,state_id,x_optmM,j,k,i_optm);
    end
end
   
figure;
col=linspecer(N*M);
tau=(0:nbin)/nbin;
list=zeros(1,2*N*M);
id=1;
for k=1:N
    for j=1:M
        list(id)=k;
        list(id+1)=j-1;
        id=id+2;
plot(tau,P(j,:,k),'-o','MarkerFaceColor',col((k-1)*N+j,:),'MarkerEdgeColor',col((k-1)*N+j,:),'MarkerSize',4);
hold on;
    end
end
hold off;
set(gca,'fontsize',14);
box on;
set(gcf,'color','w');
xlim([-0.05, 1]);
ylim([0 1.05]);
legend(sprintfc('Q%d, P%d',list));

figure;
col=linspecer(N*M);
list=zeros(1,2*N*M);
id=1;
for k=1:N
    for j=1:M
        list(id)=k;
        list(id+1)=j-1;
        id=id+2;
semilogy(tau,P(j,:,k),'-o','MarkerFaceColor',col((k-1)*N+j,:),'MarkerEdgeColor',col((k-1)*N+j,:),'MarkerSize',4);
hold on;
    end
end
hold off;
set(gca,'fontsize',14);
box on;
set(gcf,'color','w');
xlim([-0.1, 1]);
legend(sprintfc('Q%d, P%d',list));
figure;
scatter3(r(:,1),r(:,2),r(:,3),5,'filled');
box on;
set(gcf,'color','w');
set(gca,'fontsize',14);
ax.BoxStyle = 'full';