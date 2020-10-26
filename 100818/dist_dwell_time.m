load('E:/dwell/0/20.mat')
dwell_list=reshape(dwell(:,1:20),1,[]);

figure
% plot(dwell')
% histogram(dwell_list,30,'Normalization','probability')
% hold on
dwell_list(dwell_list>40)=[];
histfit(dwell_list,300,'kernel')
% [f,x]=hist(dwell_list,300);
% plot((x),(f/trapz(x,f)),'.')
% xlabel('dwell time')
% ylabel('Distribution')

%%
load('E:/dwell/0/23.mat')
con=condu(:,1);
figure
histogram(log(con),100,'Normalization','probability')
xlabel('log(T)')
ylabel('Distribution')

%%  only 25.mat
load('E:/dwell/0/25.mat')
tmp=[];
for i=1:size(dwell,1)
a=squeeze(dwell(1,:,:));
tmp=[tmp,a];
end
figure
histogram(tmp(1,:),30,'Normalization','probability')
figure
histogram(tmp(2,:),30,'Normalization','probability')
%%
%%after and include 26.mat
% load('E:/dwell/0/28.mat');
% a1=dwell;
% load('E:/dwell/0/29.mat');
% a2=dwell;
% load('E:/dwell/0/30.mat');
% a3=dwell;
% dwell=[a1; a2; a3];
t0=60/pi;
load('E:/dwell/0/27.mat');
dw=reshape(dwell(:,1:20),[],1);
dw_filter=real(dw).*(abs(imag(dw))<1e-3);
t=reshape(dwell(:,21:40),[],1);
% temp=[dw_filter t];
% save('E:/dwell/0/27_0.mat','temp')

dw_filter(dw_filter>160 )=[];
dw_filter(dw_filter<0.1 )=[];
dw_filter=dw_filter/t0;
figure
% h=histfit(dw_filter,400,'kernel');
% delete(h(2))
histogram(dw_filter,400,'Normalization','pdf','FaceAlpha',1,'EdgeColor','None')
xlabel('Time delay t (t_{0})')
ylabel('Probability Density')
set(gca,'FontSize',15)
xlim([0,6])

figure
[f,x]=hist(dw_filter,400);
plot(log(x),movmean(log(f/trapz(x,f)),5),'o')

tmp_x=log(x);
tmp_y=movmean(log(f/trapz(x,f)),5);
P = polyfit(tmp_x(40:200),tmp_y(40:200),1);
yfit = P(1)*tmp_x(40:end)+P(2);
hold on;
plot(tmp_x(40:end),yfit,'r-.')

xlabel('ln(t)')
ylabel('ln(P(t))')
legend('','P(t)~1/t^{2.18}')
legend boxoff   
set(gca,'FontSize',20)
% xlim([2.5,5.2])

t(t>1.2)=[];
t(t<0 )=[];
figure
% h=histfit(t,400,'kernel');
% delete(h(2))
% histogram(t,200,'Normalization','probability','FaceAlpha',1)
histogram(t,400,'Normalization','pdf','FaceAlpha',1,'EdgeColor','None')
xlabel('Transmission')
ylabel('Probability Density')
set(gca,'FontSize',15)
xlim([0.97,1.3])
yticks([0 70 140])
xticks([1 1.1 1.2 1.3])
% figure
% histogram(dwell(:,end),100,'Normalization','probability')

%%  
% distribution of conductance

load('50.mat')
% c1=condu;
% load('32.mat')
% c2=condu;
% load('33.mat')
% c3=condu;
% c=[c2(:,1:2); c3];

% load('33.mat')
c=condu;

figure

yyaxis left
histogram(c(:,2),100,'Normalization','pdf','FaceAlpha',1,'EdgeColor','None')
yticks([0 10 20])
% xlabel('T_{\downarrow\uparrow}','FontWeight','bold')
ylabel('Probability Density of T_{\downarrow\uparrow}')
set(gca,'FontSize',15)
yyaxis right
% hold on
% [f,x]=hist(c(:,1),100);
% bar((x),(f/trapz(x,f)))
% h=histfit(c(:,1),100,'kernel');
% histogram(c(:,1),100,'normalization','pdf')
histogram(c(:,1),100,'Normalization','pdf','FaceAlpha',1,'EdgeColor','None')
yticks([0 4 8 12])
xlabel('T_{\downarrow\uparrow} or T_{\uparrow\uparrow}')
ylabel('Probability Density of T_{\uparrow\uparrow}')
set(gca,'FontSize',15)
% alpha(.5)
% delete(h(1))
% alpha(.5)
%% after 55
load('0/88.mat')
c=condu;
figure
subplot(1,2,1)
plot(c(:,1),'o')
hold on
plot(c(:,2),'s')
ylabel('T')
xlabel('Sample Index')
yticks([0 0.25 0.5 0.75 1])
set(gca,'FontSize',10)
subplot(1,2,2)
hold on
plot(real(c(:,3)).*(abs(imag(c(:,3)))<1e-4),'o')
plot(real(c(:,4)).*(abs(imag(c(:,4)))<1e-7),'s')
ylabel('Time delay')
xlabel('Sample Index')
box on
set(gca,'FontSize',10)
% figure
% histogram(real(c(:,4)).*(abs(imag(c(:,4)))<1e-4))
%%   dos
% dos=csvread('E:\dwell\6.csv',1,1); 
% figure
% histogram(dos,20,'Normalization','pdf')
% title(['dos=',num2str(mean(dos))])

figure(10)
hold on
for m=[7,11,10,9]
    dos=csvread(['E:\dwell\',num2str(m),'.csv'],1,1); 
    dos(dos>250)=[];
    [f,x]=hist(dos,30);
    plot(x*2*sqrt(3)/4,f/trapz(x,f))    %f/trapz(x,f)
end
xlabel('DOS')
ylabel('Probability Density')
set(gca,'FontSize',15)
box on
yticks([0 0.1 0.2 0.3])
xlim([0,200])
% figure
% plot(dos)
%% lnI  after dwell/19
figure(1)
hold on
for m=601 %[60,64] %[19,20,21,22] %[60,61] %[40,38,39] %[44,43,41,40]%[40,41] % [23,29]  %[19,20,21,22]  %[25,26,27] % 28 %
csvread(['E:\dwell\',num2str(m),'.csv'],1,1); 


x=ans(:,1); y1=movmean(mean(ans(:,2:1:1000),2),8);
plot(x(1:1:end),y1(1:1:end),'.')
%     P = polyfit(x(50:350),y1(50:350),1);
%     yfit = P(1)*x+P(2);
%     hold on;
%     plot(x,yfit,'r-.','linewidth',2);
% plot(ans(:,1),movmean(mean(ans(:,3:2:end),2),2))

end
% legend('lead=0','lead=1','lead=5')
% legend('\alpha=0','\alpha=0.3','\alpha=0.5','\alpha=1')
% legend boxoff
xlabel('Position')
ylabel('<lnW(x)>')
set(gca,'FontSize',20)
box on
% xlim([-105,105])
%%
% figure(3)
% hold on

for m=68 %[60,64] %[19,20,21,22] %[40,38,39] %[44,43,41,40]%[40,41] % [23,29]  %[19,20,21,22]  %[25,26,27] % 28 %
csvread(['E:\dwell\',num2str(m),'.csv'],1,1); 
x=ans(:,1);
y=ans(:,2:end);

figure(16)
hold on

y1=movmean(mean(y(:,1:2:end),2),4);

plot(x(1:1:end),y1(1:1:end),'.')

end

%%
clear
y=[];
for t=[400,600]
csvread(['E:\dwell\',num2str(t),'.csv'],1,1);
y1=ans(:,2:1:end);
y=[y,y1];
end
x=ans(:,1);
%%
x_m=x(1:2:end);
y_m=movmean(mean(y(:,1:1:10000),2),1);
y_m_ave=mean(reshape(y_m,2,[]),1);
figure(11)
hold on
plot(x_m(1:1:end),movmean(y_m_ave(1:1:end),6),'.')
%%  70-72
g=[];
trans=[];
k=1;
for i=0:499
load(['E:\dwell\72\',num2str(i),'.mat'])
g(k)=gf(1,1);
trans(k)=t;
k=k+1;
end

figure
plot(trans)
title('T')

figure
hold on
plot(abs(g).^2,'.')
plot(abs(g).^2)