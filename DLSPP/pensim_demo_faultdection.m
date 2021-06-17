close all
clc
clear
%% load data 
load('./data/traindata.mat')
% test= load('./data/Fault3.txt'); %load fault3�����ʲ����ٶ���0.1%��б�����ӣ� ,��ʼ����������Ϊ107L
test=load('./data/Fault4.txt');% fault 4 ���������ʽ�Ծ��������10%
test = test(:,[2:5,7,10:14]);
test = reshape(test,1,size(test,1),size(test,2));
%% raw data plot 
figure
for i=1:size(train,3)
    subplot(4,3,i)
    plot(train(:,:,i)','k')
    hold on 
    plot(test(:,:,i)','r')   
end
%% Data Preprocessing: noramlization
% noramlization of train data
[final,avg,var] = zscore_batch(permute(train,[1 3 2]),2);
% noramlization of test data
test = permute(test,[1 3 2 ]);
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
%% Data Preprocessing: construct augmented data
% ѵ�����ݹ����������
final=permute(final,[3,2,1]);
temp=[];
lag=2;%ϵͳ����
for i=1:size(train,1)
    temp(:,:,i)=constructAM(final(:,:,i),lag);
end
final=permute(temp,[3,2,1]);
clear temp
% �������ݹ����������
test=permute(test,[3,2,1]);
temp=[];
lag=2;%ϵͳ����
for i=1:size(test,3)
    temp(:,:,i)=constructAM(test(:,:,i),lag);
end
test=permute(temp,[3,2,1]);
clear temp
%% ��ͼ
dd=permute(final,[1,3,2]);
cc = permute(test,[1,3,2]);
figure
for i=1:size(train,3)
    subplot(4,3,i)
    plot(dd(:,:,i)')
    hold on 
    plot(cc(:,:,i)','r')   
end
%% choose principal component
temp=permute(final,[2,1,3]);
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
[segment,~,wc] = pcaseg_batch_bu(final,10,q,1);% lag
%% �и��Ч��ͼ
figure
set(gcf,'color','w')
% 1�� ͨ���� 2�����蹦�� 3�����ʲ����ٶ� 4�����ʲ����¶�  6:DO conc �ܽ���Ũ�� 9�� �������ݻ� 10��������̼Ũ�� 11��PHֵ 12�����͹��¶� 13���������� 
legendarray={'Aeration rate','Agitator power','Substrate feed flow rate','Substrate feed temperature','DO conc','Culture volume','CO2 conc','pH','Generated heat','Fermenter temperature'};
for i=1:size(train,3)
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
%% ��׶�PCA��ģ
lx=[segment(:).lx];
rx=[segment(:).rx];
[P,FAI,num_pc,lamda,T2UCL1,QUCL1,kesi]=multiphase_pca_modelling(final,0.90,0.95,lx,rx);%95%������
[~,~,~,~,T2UCL2,QUCL2,~]=multiphase_pca_modelling(final,0.90,0.99,lx,rx);%99%������
%% �����޼���
%95%������
T2UCL_95=[];
for i=1:length(lx)
    T2UCL_95=[T2UCL_95 repmat(T2UCL1(i),1,(rx(i)-lx(i)+1))];
end
QUCL_95=[];
for i=1:length(lx)
    QUCL_95=[QUCL_95 repmat(QUCL1(i),1,(rx(i)-lx(i)+1))];
end
%99%������
T2UCL_99=[];
for i=1:length(lx)
    T2UCL_99=[T2UCL_99 repmat(T2UCL2(i),1,(rx(i)-lx(i)+1))];
end
QUCL_99=[];
for i=1:length(lx)
    QUCL_99=[QUCL_99 repmat(QUCL2(i),1,(rx(i)-lx(i)+1))];
end
%% ����
[T2,Q, fai1]=multiphase_pca_test(test,P,num_pc,lamda,T2UCL_99,QUCL_99,lx,rx);
[~,~, fai2]=multiphase_pca_test(test,P,num_pc,lamda,T2UCL_95,QUCL_95,lx,rx);
%%  plot of onlining monitoring
figure
subplot(2,1,1)
semilogy(T2','ko--','MarkerSize',5)
hold on
semilogy(T2UCL_99,'r--','linewidth',2.5)
hold on
semilogy(T2UCL_95,'k--','linewidth',2.5)

% title('BUPP Monitor Performance')
% title('DLSSPP���Ч��')
% semilogy(T2','k-','linewidth',2)
% legend('99% upper limit','95% upper limit','Location','northwest')
legend('Hotelling''s T^2','99% upper limit','95% upper limit','Location','best')

xlabel('ʱ��(h)');
ylabel('T^2');
subplot(2,1,2)
semilogy(Q','ko--','MarkerSize',5)
% semilogy(Q','k-','linewidth',2)
hold on
semilogy(QUCL_95,'k--','linewidth',2.5)
hold on
semilogy(QUCL_99','r-','linewidth',2.5)
legend('Hotelling''s T^2','99% upper limit','95% upper limit','Location','best')
set(gcf,'color','w')
xlabel('ʱ��(h)');
ylabel('Q');
set(gcf,'unit','centimeters','position',[8 5 18 12])
% legend('99% upper limit','95% upper limit','Location','northeast')