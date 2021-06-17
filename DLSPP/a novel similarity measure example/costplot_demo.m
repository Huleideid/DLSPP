clear
clc
close all
%% 
n=500;          %产生随机数的个数 
mu=[0 0];%均值
Sigma1=[1,1;1,1]; %第一组数据的协方差矩阵
Sigma2=Sigma1; %第二组数据的协方差矩阵
thea=pi/6;%每次旋转的角度
rotate_matxix=[cos(thea),-sin(thea) ;sin(thea),cos(thea)]; %旋转矩阵
level = 3;
%% 
figure
set(gcf,'color','w')
data1=mvnrnd(mu,Sigma1,n); 
data2 = data1;
for i=1:4
    subplot(2,2,i)
    scatter(data2(:,1),data2(:,2),'r'); 
    hold on
    scatter(data1(:,1),data1(:,2),'g'); 
    hold on
%     myelipsnorm(mu,Sigma1,level,'k--',2.5)
    myelipsnorm(mu,Sigma1,level,'k',2.5)
    hold on
%     myelipsnorm(mu,Sigma2,level,'k--',2.5)
    myelipsnorm(mu,Sigma2,level,'k',2.5)
    hold on
    data=[data1;data2];
    Sigma =cov(data);
    [~,ei,~]=svd(Sigma);
    elipsnorm5(mu,Sigma,level,thea/2*(i-1),1) %旋转
    myelipsnorm(mu,Sigma,level,'k-',1)  % plot control ellipse 
    elipsnorm2(mu,Sigma,level,1)% plot pc direction
    data2 = data2 * rotate_matxix;
    Sigma2=rotate_matxix'*Sigma2*rotate_matxix;
    xlabel('Var1')
    ylabel('Var2')
    axis([-4 4 -4 4])
end
