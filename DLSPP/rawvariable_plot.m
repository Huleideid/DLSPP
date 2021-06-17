%% 
load('./data/traindata.mat')
% test= load('./data/Fault3.txt'); %load fault3（基质补料速度以0.1%的斜率增加） ,初始培养基容量为107L
Fault4=load('./data/Fault4.txt');% fault 4 搅拌器速率阶跃性增长了10%
Fault4 = Fault4(:,[2:5,7,10:14]);
Fault4 = reshape(Fault4,1,size(Fault4,1),size(Fault4,2));
Fault3=load('./data/Fault3.txt');% fault 4 搅拌器速率阶跃性增长了10%
Fault3 = Fault3(:,[2:5,7,10:14]);
Fault3 = reshape(Fault3,1,size(Fault3,1),size(Fault3,2));
%% 
legendarray={'通风率(L/h)','搅拌功率(W)','基质补料速度(L/h)','基质补料温度(K)','溶解氧浓度(%)','培养基容积(L)','二氧化碳浓度(mmole/L)','pH值','发酵罐温度(K)','产生热量(kcal/h)'};
figure
for i=1:size(train,3)
    subplot(4,3,i)
    plot(train(:,:,i)','k','linewidth',2)
    hold on 
    plot(Fault3(:,:,i)','g','linewidth',2)   
    hold on 
    plot(Fault4(:,:,i)','r')   
    legend('正常数据','故障3','故障4')
    title(legendarray(i),'FontSize',10)
end
%% 
legendarray={'通风率(L/h)','搅拌功率(W)','基质补料速度(L/h)','基质补料温度(K)','溶解氧浓度(%)','培养基容积(L)','二氧化碳浓度(mmole/L)','pH值','发酵罐温度(K)','产生热量(kcal/h)'};
figure
for i=1:size(train,3)
    subplot(4,3,i)
    plot(train(:,:,i)','linewidth',2)
    title(legendarray(i),'FontSize',10)
end
set(gcf,'color','w')