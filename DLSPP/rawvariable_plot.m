%% 
load('./data/traindata.mat')
% test= load('./data/Fault3.txt'); %load fault3�����ʲ����ٶ���0.1%��б�����ӣ� ,��ʼ����������Ϊ107L
Fault4=load('./data/Fault4.txt');% fault 4 ���������ʽ�Ծ��������10%
Fault4 = Fault4(:,[2:5,7,10:14]);
Fault4 = reshape(Fault4,1,size(Fault4,1),size(Fault4,2));
Fault3=load('./data/Fault3.txt');% fault 4 ���������ʽ�Ծ��������10%
Fault3 = Fault3(:,[2:5,7,10:14]);
Fault3 = reshape(Fault3,1,size(Fault3,1),size(Fault3,2));
%% 
legendarray={'ͨ����(L/h)','���蹦��(W)','���ʲ����ٶ�(L/h)','���ʲ����¶�(K)','�ܽ���Ũ��(%)','�������ݻ�(L)','������̼Ũ��(mmole/L)','pHֵ','���͹��¶�(K)','��������(kcal/h)'};
figure
for i=1:size(train,3)
    subplot(4,3,i)
    plot(train(:,:,i)','k','linewidth',2)
    hold on 
    plot(Fault3(:,:,i)','g','linewidth',2)   
    hold on 
    plot(Fault4(:,:,i)','r')   
    legend('��������','����3','����4')
    title(legendarray(i),'FontSize',10)
end
%% 
legendarray={'ͨ����(L/h)','���蹦��(W)','���ʲ����ٶ�(L/h)','���ʲ����¶�(K)','�ܽ���Ũ��(%)','�������ݻ�(L)','������̼Ũ��(mmole/L)','pHֵ','���͹��¶�(K)','��������(kcal/h)'};
figure
for i=1:size(train,3)
    subplot(4,3,i)
    plot(train(:,:,i)','linewidth',2)
    title(legendarray(i),'FontSize',10)
end
set(gcf,'color','w')