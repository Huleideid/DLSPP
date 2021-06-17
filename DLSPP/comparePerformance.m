%% load data
load aug_traindtwc.mat
load aug_trainupwc.mat
load aug_trainswwc.mat 

%% plot
figure
set(gcf,'color','w')
plot(aug_traindtwc(:,1:20),'g','marker','o','LineWidth',2)
grid on
hold on
xlabel('阶段数')
ylabel('基于Q统计量的全局代价函数')
plot((aug_trainupwc'),'r','marker','s','LineWidth',2)
hold on 
phase=[1 phase];
sw_trainwc=[aug_traindtwc(1) sw_trainwc];
plot(phase,sw_trainwc,'b-','marker','d','LineWidth',2)
legend('DLSPP','MPPCA','SSSP')