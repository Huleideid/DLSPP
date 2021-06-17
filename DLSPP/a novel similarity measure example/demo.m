clear
clc
close all
%% 
n=500;          %����������ĸ��� 
mu=[0 0];%��ֵ
Sigma1=[0.5,0;0,sqrt(3.5)]; %��һ�����ݵ�Э�������
Sigma2=[sqrt(1.5),0;0,sqrt(2.5)]; %�ڶ������ݵ�Э�������
% Sigma2=[sqrt(3),0;0,1]; %�ڶ������ݵ�Э�������
level = 3;
%% 

data1 = mvnrnd(mu,Sigma1,n); 
data2 = mvnrnd(mu,Sigma2,n);

figure
set(gcf,'color','w')
scatter(data2(:,1),data2(:,2),'r','filled'); 
hold on
scatter(data1(:,1),data1(:,2),'g','filled'); 
hold on
myelipsnorm(mu,Sigma1,level,'--k',2.5)
hold on
myelipsnorm(mu,Sigma2,level,'--k',2.5)
xlabel('Var1')
ylabel('Var2')
axis([-4 4 -5 5])
axis equal
grid on