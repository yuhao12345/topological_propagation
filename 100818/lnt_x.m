% run lnI_position.R first

% ln_t=csvread('E:\yuhao\lnt65.csv',1,1);  %14 15 13 16 17 18 25
% ln_t2=csvread('E:\yuhao\lnt63.csv',1,1);
% ln_t=[ln_t1,ln_t2(:,2:end)];
ln_t=csvread('E:\dwell\40_0.csv',1,1);

figure(2)
hold on
lnt_ave=mean(ln_t(:,2:1:end),2);
y1=movmean(lnt_ave,5);
x1=ln_t(:,1);
% x=x1(1:5:end); y=y1(1:5:end);
plot(x1,y1,'o')
% P2 = polyfit(x1(round(length(x1)*0.3):round(length(x1)*0.75)),y1(round(length(x1)*0.3):round(length(x1)*0.75)),1);
% yfit = P2(1)*x+P2(2);
% plot(x,yfit,'-.');
% xlabel('position','FontSize', 15)
% ylabel('<lnI(x)>','FontSize', 15)
% % xlim([-500,500])
% ylim([-14,-8])
% % title('finite TI')
% loc(7)=-1/P2(1);
% 
% figure(2)
% x3=200*(1:7);
% plot(x3,loc)
% xlabel('system length')
% ylabel('effective localization length')
% p3=polyfit(x3,loc,2);
% yfit = p3(1)*x3.^2+p3(2)*x3+p3(3);
% hold on
% plot(x3,yfit,'-.');

% figure
% hold on
% x=ln_t(:,1)-min(ln_t(:,1)); y=mean(ln_t(:,2:end),2);
% plot(x,y,'.')
% xlabel('x')
% ylabel('<W(x)>')

% P2 = polyfit(x(1:600),y(1:600),1);
% yfit = P2(1)*x+P2(2);
% plot(x,yfit,'-.');
% 
% loc2=-1/P2(1);
% legend(['1/slope=',num2str(loc)])
% legend boxoff

% g=[];
% k=1;
% for i=0:4999
% load(['E:\dwell\54\',num2str(i),'.mat'])
% g(k)=t;
% k=k+1;
% end
% mean(g)
% figure
% plot(g)
box on

ylabel('<W(x)>')
xlabel('x')
set(gca, 'FontSize', 20)