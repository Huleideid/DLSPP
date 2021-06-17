clear
clc
close all
%% 数据加载和标准化
load('./data/traindata.mat')
load('./data/testdata.mat')
%% 训练数据预处理
data=train;
data=permute(data,[1 3 2]);
[data,avg,var]=zscore_batch(data,2);
lag=2;%系统阶数
data=permute(data,[3,2,1]);
temp=[];
for i=1:size(train,1)
    temp(:,:,i)=constructAM(data(:,:,i),lag);
end
data=permute(temp,[3,2,1]);
clear temp
%% 测试数据预处理
test=permute(test,[1 3 2 ]);
x_test_number=size(test,1);
for i=1:size(test,3)
    for j=1:size(test,2)
        if var(i,j)>1*10^(-10)
            test(:,j,i)=(test(:,j,i)-repmat(avg(i,j),x_test_number,1))./repmat(var(i,j),x_test_number,1);
        else
            test(:,j,i)=(test(:,j,i)-repmat(avg(i,j),x_test_number,1));
        end
    end
end
test=permute(test,[3,2,1]);
temp=[];
lag=2;%系统阶数
for i=1:size(test,3)
    temp(:,:,i)=constructAM(test(:,:,i),lag);
end
test=permute(temp,[3,2,1]);
clear temp
%% 参数设置
q=1;
min_lenth=1;%最小阶段长度
alpha=0.07;
num_segments=10;%分割终止条件
%% MPPCA
t0=cputime;%计时起点
[segment,aug_trainupwc]=pcaseg_batch_td(data,num_segments,q);
t1=cputime-t0;%计时结束
%% 阶段划分示意图可视化
u1=0;
for i=1:(length(segment) )
    u1=u1+[zeros(segment(i).lx-1,1);ones(segment(end).rx-segment(i).lx+1,1)];
end 
figure
plot(u1,'LineWidth',2,'LineStyle','none','Marker','o')
set(gcf,'color','w')
xlabel('时间(h)','FontSize',10,'FontWeight','bold','FontName','Times New Roman')
ylabel('阶段','FontSize',10, 'FontWeight','bold','FontName','Times New Roman ')
box off

%% 保存结果
num_segments=20;%分割终止条件
[segment,aug_trainupwc]=pcaseg_batch_td(data,num_segments,q);
save aug_trainupwc.mat aug_trainupwc
%% 全局代价函数可视化
figure
plot(aug_trainupwc)
set(gcf,'color','white')
%% 多阶段PCA建模
lx=[segment(:).lx];
rx=[segment(:).rx];
[P,FAI,num_pc,lamda,T2UCL1,QUCL1,kesi1]=multiphase_pca_modelling(data,0.90,0.95,lx,rx);%95%控制限
[~,~,~,~,T2UCL2,QUCL2,kesi2]=multiphase_pca_modelling(data,0.90,0.99,lx,rx);%99%控制限
% 控制限计算
T2UCL_95=[];%95%置信限
for i=1:length(lx)
    T2UCL_95=[T2UCL_95 repmat(T2UCL1(i),1,(rx(i)-lx(i)+1))];
end
QUCL_95=[];
for i=1:length(lx)
    QUCL_95=[QUCL_95 repmat(QUCL1(i),1,(rx(i)-lx(i)+1))];
end
kesi_95=[];
for i=1:length(lx)
    kesi_95=[kesi_95 repmat(kesi1(i),1,(rx(i)-lx(i)+1))];
end
T2UCL_99=[];%99%置信限
for i=1:length(lx)
    T2UCL_99=[T2UCL_99 repmat(T2UCL2(i),1,(rx(i)-lx(i)+1))];
end
QUCL_99=[];
for i=1:length(lx)
    QUCL_99=[QUCL_99 repmat(QUCL2(i),1,(rx(i)-lx(i)+1))];
end
kesi_99=[];
for i=1:length(lx)
    kesi_99=[kesi_99 repmat(kesi2(i),1,(rx(i)-lx(i)+1))];
end
%% 测试
[T2,Q, fai1]=multiphase_pca_test(test,P,num_pc,lamda,T2UCL_99,QUCL_99,lx,rx);
[~,~, fai2]=multiphase_pca_test(test,P,num_pc,lamda,T2UCL_95,QUCL_95,lx,rx);
%% 测试数据监控可视化
% figure
% subplot(2,1,1)
% semilogy(T2UCL_99,'r--','linewidth',2.5)
% hold on
% semilogy(T2UCL_95,'k--','linewidth',2.5)
% hold on
% semilogy(T2','ko--','MarkerSize',5)
% title('BUPP Monitor Performance')
% legend('99% upper limit','95% upper limit')
% xlabel('时间');
% ylabel('T2');
% subplot(2,1,2)
% semilogy(Q','ko--','MarkerSize',5)
% hold on
% semilogy(QUCL_95,'k--','linewidth',2.5)
% hold on
% semilogy(QUCL_99','r-','linewidth',2.5)
% set(gcf,'color','w')
% xlabel('时间');
% ylabel('SPE');
% set(gcf,'unit','centimeters','position',[8 5 18 12])
%% 误报率计算
FA_Q=zeros(2,size(Q,2));
for i=1:size(Q,2)
    FA_Q(1,i)=FA_Q(1,i)+sum(Q(:,i)>QUCL_99(i))/size(Q,1);
    FA_Q(2,i)=FA_Q(2,i)+sum(Q(:,i)>QUCL_95(i))/size(Q,1);
end
FA_Q=mean(FA_Q');
FA_T2=zeros(2,size(T2,2));
for i=1:size(Q,2)
    FA_T2(1,i)=FA_T2(1,i)+sum(T2(:,i)>T2UCL_99(i))/size(Q,1);
    FA_T2(2,i)=FA_T2(2,i)+sum(T2(:,i)>T2UCL_95(i))/size(Q,1);
end
FA_T2=mean(FA_T2'); 
FA_fai=zeros(2,size(fai1,2));
for i=1:size(Q,2)
    FA_fai(1,i)=FA_fai(1,i)+sum(fai1(:,i)>kesi_99(i))/size(Q,1);
    FA_fai(2,i)=FA_fai(2,i)+sum(fai2(:,i)>kesi_95(i))/size(Q,1);
end
FA_fai=mean(FA_fai');
FA_combine =zeros(2,size(fai1,2));
for i=1:size(Q,2)
    FA_combine(1,i)=FA_combine(1,i)+sum((Q(:,i)>QUCL_99(i))|(T2(:,i)>T2UCL_99(i)))/size(Q,1);
    FA_combine(2,i)=FA_combine(2,i)+sum((Q(:,i)>QUCL_95(i))|(T2(:,i)>T2UCL_95(i)))/size(Q,1);
end
FA_combine=mean(FA_combine'); 
FA_T2
FA_Q
FA_fai
FA_combine

%% 不控制
% FA_T2 =
% 
%     0.0282    0.0559
% 
% 
% FA_Q =
% 
%     0.0397    0.0575
% 
% 
% FA_fai =
% 
%     0.0422    0.0688
% FA_combine =
% 
%     0.0573    0.0920
%% 主元数为2时
% FA_T2 =
% 
%     0.0197    0.0400
% 
% 
% FA_Q =
% 
%     0.0370    0.0620
% 
% 
% FA_fai =
% 
%     0.0405    0.0678
% 
% 
% FA_combine =
% 
%     0.0501    0.0866
%% 主元数为1时
% 
% FA_T2 =
% 
%          0    0.0368
% 
% 
% FA_Q =
% 
%     0.0249    0.0557
% 
% 
% FA_fai =
% 
%     0.0327    0.0607
% 
% 
% FA_combine =
% 
%     0.0249    0.0753