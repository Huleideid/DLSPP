clear
clc
close all
%% 
n=500;          %产生随机数的个数 
mu=[0 0]; 
Sigma1=[1,1;1,1]; 
Sigma2=Sigma1; 
thea=pi/16;
rotate_matxix=[cos(thea),-sin(thea) ;sin(thea),cos(thea)];
level=2;
%利用getframe函数捕获当前plot画面，写入成gif文件
[A,map] = rgb2ind(frame2im(getframe),256);
imwrite(A,map,'cost_function.gif','LoopCount',65535,'DelayTime',0.05);
for i=1:33
    r1=mvnrnd(mu,Sigma1,n);
    r2=mvnrnd(mu,Sigma2,n); 
    scatter(r1(:,1),r1(:,2),'b'); 
    hold on
    scatter(r2(:,1),r2(:,2),'r'); 
    set(gcf,'color','w')
    hold on
    myelipsnorm(mu,Sigma1,level,'k--',2.5)
    hold on
    myelipsnorm(mu,Sigma1,level,'k--',2.5)
    hold on
    r=[r1;r2];
    sum(sum(r.^2))
    C =cov(r);
    [~,ei,~]=svd(C)
    myelipsnorm(mu,C,level,'k--',1) 
    elipsnorm2(mu,C,level,1)
    Sigma2=rotate_matxix'*Sigma2*rotate_matxix;
    xlabel('Var1')
    ylabel('Var2')
    axis([-4 4 -4 4])

    %利用getframe函数捕获当前plot画面，写入成gif文件
    [A,map] = rgb2ind(frame2im(getframe),256);
    imwrite(A,map,'cost_function.gif','LoopCount',65535,'DelayTime',0.05);
%     pause('on')
%     pause(0.01)
    hold off
end
