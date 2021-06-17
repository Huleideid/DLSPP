close all
clc
clear
DeepPink=[255 20 147]/255;
Green=[0 255  0]/255;
%% 加载数据
load('./data/traindata.mat') 
train = permute(train,[1 3 2]); 
[final,avg,var]=zscore_batch(train,2);
%% 互相关图: 选择滞后阶数
figure
set(gcf,'color','w')
lag=2;
for h=0:lag
    temp=permute(final,[3,2,1]);
    temp=constructAM(final(:,:,1),lag);   
    [U,S,V] = svd(cov(temp));
    T=temp*V;
    figure
    set(gcf,'color','w')
    for i=1:2*(h+1)
        for j=1:i
            subplot(2*(h+1),2*(h+1),(i-1)*2*(h+1)+j)
            crosscorr(T(:,i),T(:,j),[],30)
            title('')
            ylabel('')
            xlabel('')
        end
    end
end