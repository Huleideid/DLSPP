clear
clc
close all
%% 加载数据
load('./data/traindata.mat')% 加载训练数据
load('./data/testdata.mat') %加载测试数据
%% 训练数据标准化
data=train;
data=permute(data,[1 3 2]);
[data,avg,var]=zscore_batch(data,2);
%% 添加滞后向量
data=permute(data,[3,2,1]);
temp=[];
lag=2;%动态阶数
for i=1:size(data,3)
    temp(:,:,i)=constructAM(data(:,:,i),lag);
end
data=permute(temp,[3,2,1]);
%% 测试数据预处理
% 标准化
test=permute(test,[1 3 2 ]);
x_test_number=size(test,1);
for i=1:size(test,3)%
    for j=1:size(test,2)
        if var(i,j)>1*10^(-10) % small trick: if var is very samll, data will only substaract mean
            test(:,j,i)=(test(:,j,i)-repmat(avg(i,j),x_test_number,1))./repmat(var(i,j),x_test_number,1);
        else
            test(:,j,i)=(test(:,j,i)-repmat(avg(i,j),x_test_number,1));
        end
    end
end

%添加滞后向量
test=permute(test,[3,2,1]);
temp=[];
lag=2;%系统阶数
for i=1:size(test,3)
    temp(:,:,i)=constructAM(test(:,:,i),lag);
end
test=permute(temp,[3,2,1]);
clear temp
%% 初始化参数
q=1; %主元数
min_lenth=3; %
sw_trainwc=[];
sw_testwc=[];
phase=[];
%% phase numer 2
[train_segment, wc] =pcaseg_batch_sw(data,min_lenth,q,4);
phase=[phase length(train_segment)];
% wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
wc = pcaseg_batch_sw_evaluate(test,q,train_segment);
sw_testwc=[sw_testwc wc];
%% phase numer 3
min_lenth=4; %
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.5);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
wc = pcaseg_batch_sw_evaluate(test,q,train_segment);
sw_testwc=[sw_testwc wc];
%% phase numer 4
min_lenth=5; 
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.1);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
wc = pcaseg_batch_sw_evaluate(test,q,train_segment);
sw_testwc=[sw_testwc wc];
%% phase numer 5 
min_lenth=3; 
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.39);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
wc = pcaseg_batch_sw_evaluate(test,q,train_segment);
sw_testwc=[sw_testwc wc];
%% phase numer 7
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.3);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
wc = pcaseg_batch_sw_evaluate(test,q,train_segment);
sw_testwc=[sw_testwc wc];
%% phase numer 8
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.1);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
wc = pcaseg_batch_sw_evaluate(test,q,train_segment);
sw_testwc=[sw_testwc wc];
%% phase numer 9
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.19);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
wc = pcaseg_batch_sw_evaluate(test,q,train_segment);
sw_testwc=[sw_testwc wc];
%% phase numer 10
t0=cputime;%计时起点
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.15);
t1=cputime-t0;%计时结束
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
%% phase numer 11
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.17);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
%% phase numer 12
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.02);
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
%% phase numer 13
t0=cputime;%计时起点
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.05);
t1=cputime-t0;%计时结束
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];
%% phase numer 16
min_lenth=2;
t0=cputime;%计时起点
train_segment=pcaseg_batch_sw(data,min_lenth,q,1.01);
t1=cputime-t0;%计时结束
phase=[phase length(train_segment)];
wc = pcaseg_batch_sw_evaluate(data,q,train_segment);
sw_trainwc=[sw_trainwc wc];

%% phase numer 10
min_lenth=4;
t0=cputime;%计时起点
train_segment=pcaseg_batch_sw(data,min_lenth, q, 1.05);
t1=cputime-t0;%计时结束
%% 阶段划分示意图
u1=0;
for i=1:(length(train_segment) )
    u1=u1+[zeros(train_segment(i).lx-1,1);ones(train_segment(end).rx-train_segment(i).lx+1,1)];
end 
figure(1)
plot(u1,'LineWidth',2,'LineStyle','none','Marker','o')
set(gcf,'color','w')
xlabel('时间(h)','FontSize',10,'FontWeight','bold','FontName','Times New Roman')
ylabel('阶段','FontSize',10, 'FontWeight','bold','FontName','Times New Roman ')
box off
%% 存储数据for监控
lx=[train_segment(:).lx];
rx=[train_segment(:).rx];
% [~,~,num_pc,~,~,~,~]=subpca_train(data,0.90,0.95,lx,rx);%95%控制限
% save sw.mat num_pc,lx,rx;
%% 
[P,FAI,num_pc,lamda,T2UCL1,QUCL1,kesi1]=multiphase_pca_modelling(data,0.90,0.95,lx,rx);%95%控制限
[~,~,~,~,T2UCL2,QUCL2,kesi2]=multiphase_pca_modelling(data,0.90,0.99,lx,rx);%99%控制限
%% 控制限计算
%95%置信限
T2UCL_95=[];
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
%99%置信限
T2UCL_99=[];
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
% [T2,Q]=subpca_test(test,P,num_pc,lamda,lx,rx);
[T2,Q, fai1]=multiphase_pca_test(test,P,num_pc,lamda,T2UCL_99,QUCL_99,lx,rx);
[~,~, fai2]=multiphase_pca_test(test,P,num_pc,lamda,T2UCL_95,QUCL_95,lx,rx);

%% 
figure
subplot(2,1,1)
semilogy(T2UCL_99,'r--','linewidth',2.5)
hold on
semilogy(T2UCL_95,'k--','linewidth',2.5)
hold on
semilogy(T2','ko--','MarkerSize',5)
title('SSSP Monitor Performance')
% semilogy(T2','k-','linewidth',2)
legend('99% upper limit','95% upper limit')
xlabel('时间');
ylabel('T2');
subplot(2,1,2)
% plot(Q','k')
semilogy(Q','ko--','MarkerSize',5)
% semilogy(Q','k-','linewidth',2)
hold on
semilogy(QUCL_95,'k--','linewidth',2.5)
hold on
semilogy(QUCL_99','r-','linewidth',2.5)
set(gcf,'color','w')
xlabel('时间');
ylabel('SPE');
set(gcf,'unit','centimeters','position',[8 5 18 12])
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
% 不控制主元的情况下
% FA_T2 =
% 
%     0.0322    0.0596
% 
% 
% FA_Q =
% 
%     0.0349    0.0541
% 
% 
% FA_fai =
% 
%     0.0406    0.0666
% 
% 
% FA_combine =
% 
%     0.0552    0.0924


% 主元数为2的情况下
% FA_T2 =
% 
%     0.0189    0.0432
% 
% 
% FA_Q =
% 
%     0.0344    0.0559
% 
% 
% FA_fai =
% 
%     0.0422    0.0683
% 
% 
% FA_combine =
% 
%     0.0478    0.0843
%% 
 figure
 plot(phase,sw_trainwc)
%% 
save aug_trainswwc.mat phase sw_trainwc