clear
path='C:/Users/ykang/Documents/correlated_loc_downdstair/36/';
t=[];
dt=[];
for i=0:19
    load([path,num2str(i),'.mat']);
    t(:,i+1)=mean(tran,1);
    dt(:,i+1)=std(tran,1);
end
surf(t)
plot(t)
xlabel('disorder strength')
ylabel('<G>')
figure
surf(dt)
plot(dt)
xlabel('disorder strength')
ylabel('std(G)')
