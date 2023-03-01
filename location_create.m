clear all
N=100;%节点个数
xm=100;ym=100;%监测范围
location=zeros(N,2);%位置
for i=1:N
    location(i,:)=[rand*xm rand*ym];%节点位置
end
save location