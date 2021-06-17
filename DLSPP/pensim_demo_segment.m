close all
clc
clear
DeepPink=[255 20 147]/255;
Green=[0 255  0]/255;
%% ��������
load('./data/traindata.mat')
load('./data/testdata.mat')
test = permute(test,[1 3 2 ]); 
train = permute(train,[1 3 2]); 
[final,avg,var]=zscore_batch(train,2);
%% �����ͼ,ѡ���ͺ����
% figure
% set(gcf,'color','w')
% lag=4;
% for h=0:lag
%     temp=permute(final,[3,2,1]);
%     temp=constructAM(final(:,:,1),lag);   
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
%% �����������
final=permute(final,[3,2,1]);
temp=[];
lag=2;%ϵͳ����
for i=1:size(train,1)
    temp(:,:,i)=constructAM(final(:,:,i),lag);
end
final=permute(temp,[3,2,1]);
clear temp
%% ��ͼ
% dd=permute(final,[1,3,2]);
% figure
% for i=1:size(train,3)
%     subplot(4,3,i)
%     plot(dd(:,:,i)')
% end
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
%% ����Ŀ��ָ���ĿΪ2����¼ȫ�ִ��ۺ�����׶���Ŀ�ı仯����
num=50;
[~,~,G_T2] = pcaseg_batch_bu(final,2,q,0);
[~,~,G_Q] = pcaseg_batch_bu(final,2,q,1);
data_vfold=reshape(permute(final,[2 1 3]),size(final,2),[]);
wcQ0 = pcaresid(data_vfold',q,1);
wcT0 = pcaresid(data_vfold',q,0);
G_T2=fliplr(G_T2');%ת��
G_T2=[wcT0,G_T2];
G_Q=fliplr(G_Q');%ת��
G_Q=[wcQ0,G_Q];
G_T2=G_T2(1:num);
G_Q=G_Q(1:num);
G = G_T2 + G_Q;%sum of Wc_T2 and Wc_Q;Theoretically, it should be constant
xlabel('ʱ�θ���')
%% ȫ�ִ��ۺ������ӻ�
figure
subplot(2,1,1)
set(gcf,'color','w')
phaserange=1:length(G_T2);
plot(phaserange,G_T2,'k-','Marker','>','MarkerSize',4.5,'color',Green,'LineWidth',1.5)
ylabel('G_{T^2}')
xlabel('phase number')
box off
subplot(2,1,2)
plot(phaserange,G_Q,'k-','Marker','<','MarkerSize',4.5,'color',DeepPink,'LineWidth',1.5)
% axis([phaserange 0.3 0.75]);
ylabel('G_Q')
xlabel('phase number')% xlabel('ʱ�θ���')
box off
%% ����Qȫ�ִ��ۺ����仯�����ӻ�
delta_SPE_cost=G_Q(1:end-1)-G_Q(2:end);
figure
% plot(delta_SPE_cost,'b','LineWidth',2.5)
plot(delta_SPE_cost,'b--o','LineWidth',2.5)
grid on
axis([0 50 0 0.2])
set(gcf,'color','w')
xlabel('�׶�����','FontSize',12)% xlabel('segment number')
ylabel('����Qͳ������ȫ�ִ��ۺ����仯��','FontSize',12)% ylabel('delta cost')
title('Qͳ������ȫ�ִ��ۺ����仯��','FontSize',12)
%% ������, ����T2��Qͳ������ȫ�ִ��ۺ�����Ӧ��Ϊ����
figure
plot(G,'b')
set(gcf,'color','w')
xlabel('�׶�����','FontSize',11)
ylabel('����T2��Qͳ������ȫ�ִ��ۺ�����','FontSize',11)
%% relate rate curve with the segment number
% rela_rate_T2=(G_T2(2:end)-G_T2(1:end-1))./G_T2(1:end-1);
% rela_rate_SPE=(G_Q(1:end-1)-G_Q(2:end))./G_Q(1:end-1);
% figure
% subplot(2,1,1)
% plot(rela_rate_T2)
% subplot(2,1,2)
% plot(rela_rate_SPE)
% set(gcf,'color','w')
%% DLSSPP
t0=cputime;
% [segment,~,wc] = pcaseg_batch_bu(final,4,q,1);%
[segment,~,wc] = pcaseg_batch_bu(final,10,q,1);
t1=cputime-t0;
%% �и����ԭʼ�����ϵ�Ч��ͼ
dd=permute(train,[1,3,2]);
figure
set(gcf,'color','w')
% 1�� aeration rate ͨ���� 2�����蹦�� 3�����ʲ����ٶ�4�����ʲ����¶�  6:DO conc �ܽ���Ũ�� 9�� �������ݻ� 10��������̼Ũ�� 11��PHֵ 12�����͹��¶� 13���������� 
legendarray={'Aeration rate','Agitator power','Substrate feed flow rate','Substrate feed temperature','DO conc','Culture volume','CO2 conc','pH','Generated heat','Fermenter temperature'};
legendarray={'ͨ����(L/h)','���蹦��(W)','���ʲ����ٶ�(L/h)','���ʲ����¶�(K)','�ܽ���Ũ��(%)','�������ݻ�(L)','������̼Ũ��(mmole/L)','pHֵ','���͹��¶�(K)','��������(kcal/h)'};
for i=1:size(train,2)
    subplot(4,3,i)
    plot(dd(:,:,i)')
    hold on
    b=axis;
    b=b(3:4);
    title(legendarray(i),'FontSize',10)
    xlabel('Time(h)','FontSize',10)
%     ylabel(legendarray(i),'FontSize',11)
end
%% �и����ԭʼ�����ϵ�Ч��ͼ
dd=permute(final,[1,3,2]);
figure
set(gcf,'color','w')
% 1�� aeration rate  2�����蹦�� 3�����ʲ����ٶ�4�����ʲ����¶�  6:DO conc �ܽ���Ũ�� 9�� �������ݻ� 10��������̼Ũ�� 11��PHֵ 12�����͹��¶� 13���������� 
legendarray={'Aeration rate','Agitator power','Substrate feed flow rate','Substrate feed temperature','DO conc','Culture volume','CO2 conc','pH','Generated heat','Fermenter temperature'};
for i=1:size(train,2)
    subplot(4,3,i)
    plot(dd(:,:,i)')
    hold on
    b=axis;
    b=b(3:4);
    title(legendarray(i))
    for j=1:size(segment,2)
        line([segment(j).lx  segment(j).lx], b,'color','g','linewidth',2.5);
    end
end
%% �׶λ���ʾ��ͼ
% figure
% plot(u1,'LineWidth',2,'LineStyle','none','Marker','o')
% set(gcf,'color','w')
% % xlabel('Time(h)','FontSize',11,'FontWeight','bold','FontName','Times New Roman')
% % ylabel('The number of phases','FontSize',11, 'FontWeight','bold','FontName','Times New Roman ')
% xlabel('Time(h)','FontSize',11,'FontWeight','bold','FontName','Times New Roman')
% ylabel('The number of phases','FontSize',11, 'FontWeight','bold','FontName','Times New Roman ')
% box off
u1=0;
for i=1:(length(segment) )
    u1=u1+[zeros(segment(i).lx-1,1);ones(segment(end).rx-segment(i).lx+1,1)];
end 
figure
k= 1:segment(end).rx;
[xb,yb] = stairs(k,u1); % gives us the x and y coordinates 
stairs(k,u1,'k--','LineWidth',2); % the plot
hold on
for i = 1:2:(2*400-2) % size(xb) == size(yb) = 2*n-1
    % overplot black line segments
    plot(xb(i:i+1),yb(i:i+1),'k-','LineWidth',3)
end
grid on
set(gcf,'color','w')
axis([0 400 0 10.5])
% xlabel('Time(h)','FontSize',11,'FontWeight','bold','FontName','Times New Roman')
% ylabel('The number of phases','FontSize',11, 'FontWeight','bold','FontName','Times New Roman ')
xlabel('ʱ��(h)','FontSize',11)% xlabel('ʱ��(h)','FontSize',11,'FontWeight','bold','FontName','Times New Roman')
ylabel('�׶����','FontSize',11)% ylabel('�׶����','FontSize',11, 'FontWeight','bold','FontName','Times New Roman ')
% b=axis;
% b=b(3:4);
% Tag = [15 45 292];
% for j=1:size(Tag,2)
%     line([Tag(j) Tag(j)], b,'color','g','linewidth',2,'linestyle','--');
% end

