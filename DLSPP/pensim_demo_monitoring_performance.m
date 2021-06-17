close all
clc
clear
%% load data 
load('./data/traindata.mat')
load('./data/testdata.mat')
test=permute(test,[1 3 2 ]); 
train =  permute(train,[1 3 2]);
%% �����ͼ,ѡ���ͺ����
% figure
% set(gcf,'color','w')
% lag=4;
% for h=0:lag
%     temp=permute(data,[3,2,1]);
%     temp=constructAM(data(:,:,1),lag);   
%     [U,S,V] = svd(cov(temp));
%     T=temp*V;
%     figure
%     set(gcf,'color','w')
%     for i=1:2*(h+1)
%         for j=1:i
%             subplot(2*(h+1),2*(h+1),(i-1)*2*(h+1)+j)
%             crosscorr(T(:,i),T(:,j),[],30)
%             title('')
%             ylabel('')
%             xlabel('')
%         end
%     end
% end
%% ѵ������Ԥ����
%��׼��
[data,avg,var]=zscore_batch(train,2);
%������̬�������
data=permute(data,[3,2,1]);
temp=[];
lag=2;%ϵͳ����
for i=1:size(train,1)
    temp(:,:,i)=constructAM(data(:,:,i),lag);
end
data=permute(temp,[3,2,1]);
clear temp
%% ��������Ԥ����
%��׼��
x_test_number=size(test,1);
for i=1:size(test,3)%������
    for j=1:size(test,2)
        if var(i,j)>1*10^(-10)
            test(:,j,i)=(test(:,j,i)-repmat(avg(i,j),x_test_number,1))./repmat(var(i,j),x_test_number,1);
        else
            test(:,j,i)=(test(:,j,i)-repmat(avg(i,j),x_test_number,1));
        end
    end
end
%������̬�������
test=permute(test,[3,2,1]);
temp=[];
lag=2;%ϵͳ����
for i=1:size(test,3)
    temp(:,:,i)=constructAM(test(:,:,i),lag);
end
test=permute(temp,[3,2,1]);
clear temp
%% ��ͼ
dd=permute(data,[1,3,2]);
figure
for i=1:size(train,2)
    subplot(4,3,i)
    plot(dd(:,:,i)')
end
%% ��Ԫ��ѡ��
temp=permute(data,[2,1,3]);
p=reshape(temp,size(temp,1),[]);
x=zscore(p');
rand('seed',564)
N=size(x,1);  % Number of samples
nx=size(x,2);    % Number of observable variables
inic=10;     % Initial number of segments
ncov=inic;
dumm=[1 ceil([1:ncov].*N./ncov)];   
for i=1:ncov
    eigvect=[]; eigval=[];
    [eigvect eigval]=eig(cov(x(dumm(i):dumm(i+1),:)));
    eigval=-1*(sort(-1*diag(eigval)));
    figure(3)
    subplot(2,1,1)
    plot(eigval,'-o')
    ylabel('eigenvalues')
    hold on
    title('Screeplot')
    subplot(2,1,2)
    plot(cumsum(eigval)/sum(eigval),'-o')
    xlabel('number of eigenvalues')
    ylabel('rate of cumulative sum')
    hold on
end
drawnow
disp('See the screeplot.')
q = input('Number of principal components will be ');
%% 
for flag=0:1
    [~,Tc,Wc] = pcaseg_batch_bu(data,2,q,flag,1);
    Wc=fliplr(Wc');%ת��
    Wc=Wc(1:200);
    Tc=fliplr(Tc');%ת��
    Tc=Tc(1:200);
    figure
    subplot(2,1,1)
    set(gcf,'color','w')
    plot(Wc,'k-','Marker','s')
    ylabel('Wc')
    xlabel('ʱ�θ���')
    subplot(2,1,2)
    plot(Tc,'k-','Marker','s')
    ylabel('Tc')
    xlabel('ʱ�θ���')
end

%% 
num=50;
data_vfold=reshape(permute(data,[2 1 3]),size(data,2),[]);
wcQ0 = pcaresid(data_vfold',q,1);
wcT0 = pcaresid(data_vfold',q,0);
[~,~,T2_cost] = pcaseg_batch_bu(data,2,q,0);
[~,~,SPE_cost] = pcaseg_batch_bu(data,2,q,1);
data_vfold=reshape(permute(data,[2 1 3]),size(data,2),[]);
wcQ0 = pcaresid(data_vfold',q,1);
wcT0 = pcaresid(data_vfold',q,0);
T2_cost=fliplr(T2_cost');%ת��
T2_cost=[wcT0,T2_cost];
SPE_cost=fliplr(SPE_cost');%ת��
SPE_cost=[wcQ0,SPE_cost];
T2_cost=T2_cost(1:num);
SPE_cost=SPE_cost(1:num);
%% 
aug_traindtwc=SPE_cost;
save aug_traindtwc.mat aug_traindtwc
%% 
Blue1=[0 0 255]/255;
DeepPink=[255 20 147]/255;
Orange=[255 165 0]/255;
Green=[0 255  0]/255;
Firebrick1=[255 48 48]/255;
DeepSkyBlue1=[0 191 255]/255;
MediumOrchid4=[122 55 139]/255;
MediumVioletRed=[199 21 133]/255;
figure
subplot(2,1,1)
set(gcf,'color','w')
phaserange=1:50;
plot(phaserange,T2_cost,'k-','Marker','>','MarkerSize',4.5,'color',Green,'LineWidth',1.5)
ylabel('T2 weightcost')
% xlabel('ʱ�θ���')
xlabel('phase number')
box off
subplot(2,1,2)
plot(phaserange,SPE_cost,'k-','Marker','<','MarkerSize',4.5,'color',DeepPink,'LineWidth',1.5)
% axis([phaserange 0.3 0.75]);
ylabel('SPE weightcost')
% xlabel('ʱ�θ���')
xlabel('phase number')
box off
%% relate rate
rela_rate_T2=(T2_cost(2:end)-T2_cost(1:end-1))./T2_cost(1:end-1);
rela_rate_SPE=(SPE_cost(1:end-1)-SPE_cost(2:end))./SPE_cost(1:end-1);
figure
subplot(2,1,1)
plot(rela_rate_T2)
subplot(2,1,2)
plot(rela_rate_SPE)
set(gcf,'color','w')
%% �и��
[segment,~,wc] = pcaseg_batch_bu(data,10,q,1);% lag
%% �и�����������ʾ��ͼ
figure
set(gcf,'color','w')
% 1�� aeration rate ͨ���� 2�����蹦�� 3�����ʲ����ٶ�4�����ʲ����¶�  6:DO conc �ܽ���Ũ�� 9�� �������ݻ� 10��������̼Ũ�� 11��PHֵ 12�����͹��¶� 13���������� 
legendarray={'Aeration rate','Agitator power','Substrate feed flow rate','Substrate feed temperature','DO conc','Culture volume','CO2 conc','pH','Generated heat','Fermenter temperature'};
for i=1:size(train,2)
    subplot(4,3,i)
    plot(dd(:,:,i)')
    hold on
    b=axis;
    b=b(3:4);
     title(legendarray(i))
    for j=1:size(segment,2)
        line([segment(j).lx  segment(j).lx], b,'color','g','linewidth',2);
    end
end
%% �׶λ���ʾ��ͼ
u1=0;
for i=1:(length(segment) )
    u1=u1+[zeros(segment(i).lx-1,1);ones(segment(end).rx-segment(i).lx+1,1)];
end 
figure
plot(u1,'LineWidth',2,'LineStyle','none','Marker','o')
set(gcf,'color','w')
xlabel('time(h)','FontSize',10,'FontWeight','bold','FontName','Times New Roman')
ylabel('phase','FontSize',10, 'FontWeight','bold','FontName','Times New Roman ')
box off
%% ��׶ν�ģ
lx=[segment(:).lx];
rx=[segment(:).rx];
[P,FAI,num_pc,lamda,T2UCL1,QUCL1,kesi1]=multiphase_pca_modelling(data,0.90,0.95,lx,rx);%95%������
[~,~,~,~,T2UCL2,QUCL2,kesi2]=multiphase_pca_modelling(data,0.90,0.99,lx,rx);%99%������
% �����޼���
T2UCL_95=[];%95%������
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
T2UCL_99=[];%99%������
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
%% 
u1=0;
for i=1:(length(segment) )
    u1=u1+[zeros(segment(i).lx-1,1);ones(segment(end).rx-segment(i).lx+1,1)];
end 
pc = zeros(400,1);
for i=1:length(segment)
    pc(segment(i).lx:segment(i).rx) =num_pc(i);
end 
figure
plot(pc,'LineWidth',3,'Color','k')
% ,'LineStyle','--'
axis([0 400 0 5])
set(gcf,'color','w')
xlabel('Time(h)','FontSize',10,'FontWeight','bold','FontName','Times New Roman')
ylabel('Number of principal components','FontSize',10, 'FontWeight','bold','FontName','Times New Roman ') 
%% ����
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
title('DLSSPP���Ч��')
% semilogy(T2','k-','linewidth',2)
legend('99% upper limit','95% upper limit')
xlabel('ʱ��');
ylabel('T2');
subplot(2,1,2)
semilogy(Q','ko--','MarkerSize',5)
% semilogy(Q','k-','linewidth',2)
hold on
semilogy(QUCL_95,'k--','linewidth',2.5)
hold on
semilogy(QUCL_99','r-','linewidth',2.5)
set(gcf,'color','w')
xlabel('ʱ��');
ylabel('SPE');
set(gcf,'unit','centimeters','position',[8 5 18 12])
%% ���ʼ���
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
%% ��ͳһ������Ԫ��
% FA_T2 =
% 
%     0.0280    0.0554
% 
% 
% FA_Q =
% 
%     0.0397    0.0605
% 
% 
% FA_fai =
% 
%     0.0435    0.0697
% 
% 
% FA_combine =
% 
%     0.0561    0.0926
%% ��Ԫ��Ϊ2ʱ
% FA_T2 =
% 
%     0.0175    0.0382
% 
% 
% FA_Q =
% 
%     0.0385    0.0619
% 
% 
% FA_fai =
% 
%     0.0404    0.0670
% 
% 
% FA_combine =
% 
%     0.0495    0.0853
%% ��Ԫ��Ϊ1ʱ
% FA_T2 =
% 
%          0    0.0375
% 
% 
% FA_Q =
% 
%     0.0340    0.0610
% 
% 
% FA_fai =
% 
%     0.0361    0.0622
% 
% 
% FA_combine =
% 
%     0.0340    0.0798