clear
width=49;
length=159;
shape_whole=[];
shape_tau1=[];
v_dist=[];
k=1;
m=1;
save_path = 'E:/correlation pattern/10/';
% load([save_path,'velocity.mat'])
for cishu =0:499
    load([save_path,num2str(cishu),'.mat'])
    wf_ave=mean(abs(wf).^2,1);  %diag(1./mode_v)*abs(wf).^2
    wf_ave_reshape=reshape(wf_ave,width,length);
    wf_shape=sum(wf_ave_reshape,1);
    shape_whole(k,:)=wf_shape;
   
    

    
%     if s(1)>0.98
%         f_tau1=abs(sum(diag(v(1,:))'*wf,1)).^2;
%         f_tau1_reshape=reshape(f_tau1,width,length);
%         shape_tau1(m,:)=sum(f_tau1_reshape,1);      
% 
%         v_dist(m,:)=abs(v(1,:)).^2;
%         m=m+1;
%     end
    k=k+1;
end
% x=linspace(0,1,length);
figure(6)
hold on
% title('ave whole eigenchannel')
l=1:length;
temp=mean(shape_whole,1);
plot(l,temp)

c = polyfit(l,temp,1);
l0=-c(2)/c(1);

y_est = polyval(c,[l0 l]);
hold on
plot([l0 l],y_est,'r--')
xlabel('sample length')
ylabel('<W(x)>')

legend('n=1 zb=18.45','n=1.1 zb=19.16','n=1.2 zb=22.63','n=1.5 zb=34.42')
% figure(2)
% hold on
% x=(0:length-1)/(length-1);
% plot(x,mean(shape_tau1,1),'-')
% legend('l=159','l=99','l=39')
% ylabel('average of tau=1 eigenchannel')
% xlabel('x/l')
% title('ave tau1')
% figure(3)
% hold on
% v_mean=mean(v_dist,1);
% plot(v_mean)
% title('ave v tau1')
% figure
% % 1/mean(1./mode_v)
% wf_1_reshape=reshape(wf(4,:),width,length);
% wf1_shape=sum(wf_1_reshape,1);
% plot(abs(wf1_shape).^2)

% zb=[];
l0-length