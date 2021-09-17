% loc length from transfer matrix, conductance g
clear
file_path='E:/loc scale/2/';
s_list=[];
t_list=[];
k=1;
for i=1:9
    load([file_path,num2str(i),'.mat'])
    t00=s10-s11/s01*s00;
    t01=s11/s01;
    t10=-inv(s01)*s00;
    t11=inv(s01);
    transfer=[t00 t01
        t10 t11];
    s=svd(transfer);
    s_list(k,:)=s;
    t_list(k)=t;
    k=k+1;
end
xi_inverse=mean(log(s_list(:,1))) %/double(length);