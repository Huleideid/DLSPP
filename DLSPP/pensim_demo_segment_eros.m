close all
clc
clear
%% ��������
load('./data/traindata.mat')
load('./data/testdata.mat')
test=permute(test,[1 3 2 ]); 
[final,avg,var]=zscore_bactch(permute(train,[1 3 2]),2);
%% ����Ԥ������׼����ʽ2
final=permute(final,[3,2,1]);
temp=[];
lag=0;%ϵͳ����
for i=1:size(train,1)
    temp(:,:,i)=constructAM(final(:,:,i),lag);
end
final=permute(temp,[3,2,1]);
clear temp
%% ��ͼ
dd=permute(final,[1,3,2]);
figure
for i=1:size(train,3)
    subplot(4,3,i)
    plot(dd(:,:,i)')
end
%% ��Ԫ��ѡ��
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
[segment,~,Wc] = pcaseg_batch_eros(final,2,q,flag,1);
Wc=fliplr(Wc');%ת��
Wc=Wc(1:200);
figure
set(gcf,'color','w')
plot(Wc,'k-','Marker','s')
ylabel('Wc')
xlabel('ʱ�θ���')
%% T2/SPE����
Blue1=[0 0 255]/255;
DeepPink=[255 20 147]/255;
Orange=[255 165 0]/255;
Green=[0 255  0]/255;
Firebrick1=[255 48 48]/255;
DeepSkyBlue1=[0 191 255]/255;
MediumOrchid4=[122 55 139]/255;
MediumVioletRed=[199 21 133]/255;
color=[Blue1;Green;Firebrick1;Firebrick1;MediumOrchid4];
%% 
num=50;
[~,~,T2_cost] = pcaseg_batch_eros(final,2,q,0);
[~,~,SPE_cost] = pcaseg_batch_eros(final,2,q,1);
data_vfold=reshape(permute(final,[2 1 3]),size(final,2),[]);
wcQ0 = pcaresid(data_vfold',q,1);
wcT0 = pcaresid(data_vfold',q,0);
T2_cost=fliplr(T2_cost');%ת��
T2_cost=[wcT0,T2_cost];
SPE_cost=fliplr(SPE_cost');%ת��
SPE_cost=[wcQ0,SPE_cost];
T2_cost=T2_cost(1:num);
SPE_cost=SPE_cost(1:num);
%% ����cost֮��Ӧ���ǳ���
% cost=SPE_cost+T2_cost
%% 
delta_SPE_cost=SPE_cost(1:end-1)-SPE_cost(2:end);
%% 
figure
plot(delta_SPE_cost,'b')
axis([0 50 0 0.2])
set(gcf,'color','w')
% xlabel('segment number')
% ylabel('delta cost')
% xlabel('�׶�����','FontSize',11,'FontWeight','bold','FontName','Times New Roman')
% ylabel('ȫ�ִ��ۺ����仯��','FontSize',11, 'FontWeight','bold','FontName','Times New Roman ')
xlabel('�׶�����','FontSize',11)
ylabel('����Qͳ������ȫ�ִ��ۺ����仯��','FontSize',11)
%% 
figure
plot(SPE_cost(1:end),'k-','Marker','s');
set(gcf,'color','w')
axis([0 50 3.5 5.2])
%% 
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
%% 
t0=cputime;
[segment,~,wc] = pcaseg_batch_eros(final,10,q,flag,1);
t1=cputime-t0;

%% �и��Ч��ͼ
figure
set(gcf,'color','w')
% 1�� aeration rate ͨ���� 2�����蹦�� 3�����ʲ����ٶ�4�����ʲ����¶�  6:DO conc �ܽ���Ũ�� 9�� �������ݻ� 10��������̼Ũ�� 11��PHֵ 12�����͹��¶� 13���������� 
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
% xlabel('Time(h)','FontSize',11,'FontWeight','bold','FontName','Times New Roman')
% ylabel('The number of phases','FontSize',11, 'FontWeight','bold','FontName','Times New Roman ')
xlabel('Time(h)','FontSize',11,'FontWeight','bold','FontName','Times New Roman')
ylabel('The number of phases','FontSize',11, 'FontWeight','bold','FontName','Times New Roman ')
box off
%% 
% u1=0;
% for i=1:(length(segment) )
%     u1=u1+[zeros(segment(i).lx-1,1);ones(segment(end).rx-segment(i).lx+1,1)];
% end 
% figure
% plot(u1,'LineWidth',2')
% axis([0 400 0 12])
% set(gcf,'color','w')
% xlabel('Time(h)','FontSize',10,'FontWeight','bold','FontName','Times New Roman')
% ylabel('phase','FontSize',10, 'FontWeight','bold','FontName','Times New Roman ')
%% 
figure
k= 1:400;
[xb,yb] = stairs(k,u1); % gives us the x and y coordinates 
stairs(k,u1,'k--','LineWidth',2); % the plot
hold on
for i = 1:2:(2*400-2) % size(xb) == size(yb) = 2*n-1
    % overplot black line segments
    plot(xb(i:i+1),yb(i:i+1),'k-','LineWidth',3)
end
set(gcf,'color','w')
axis([0 400 0 10.5])
xlabel('ʱ��(h)','FontSize',11)
ylabel('�׶����','FontSize',11)

