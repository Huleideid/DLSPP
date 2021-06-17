function [P,FAI,num_pc,lamda,T2UCL1,QUCL,kesi]= pca_batch(Xtrain,threhold,alpha)
    %% 无需中心化和标准化，因为经过间歇过程预处理后可以执行了该步骤
    [X_row,X_col] = size(Xtrain); %求Xtrain行、列数               
    %求协方差矩阵
    sigmaXtrain = cov(Xtrain);
    %对协方差矩阵进行特征分解，lamda为特征值构成的对角阵，T的列为单位特征向量，且与lamda中的特征值一一对应：
    [T,lamda] = eig(sigmaXtrain);                            
                                             
    %取对角元素(结果为一列向量)，即lamda值，并上下反转使其从大到小排列，主元个数初值为1，若累计贡献率小于90%则增加主元个数
    D = flipud(diag(lamda));                                                                           
    num_pc = 1;                                         
    while sum(D(1:num_pc))/sum(D) < threhold  
        num_pc = num_pc +1;
    end   
    %取与lamda相对应的特征向量
    P = T(:,X_col-num_pc+1:X_col); 
   %% 计算控制限
    %求置信度为99%、95%时的T2统计控制限                       
    T2UCL1=num_pc*(X_row-1)*(X_row+1)*finv(alpha,num_pc,X_row - num_pc)/(X_row*(X_row - num_pc));
    %置信度为95%的Q统计控制限
        for i = 1:3
            theta(i) = sum((D(num_pc+1:X_col)).^i);
        end
    if (max(theta)>1*10^-5)       
%         g = 1/(X_row-1)*sum( D(num_pc+1:X_col).^2)/sum(D(num_pc+1:X_col));
%         h = (sum(D(num_pc+1:X_col).^2))/sum(D(num_pc+1:X_col).^2);
%         QUCL = g*chi2inv(alpha, h); 
        h0 = 1 - 2*theta(1)*theta(3)/(3*theta(1)^2);
        ca = norminv(alpha,0,1);
        QUCL = theta(1)*(h0*ca*sqrt(2*theta(2))/theta(1) + 1 + theta(2)*h0*(h0 - 1)/theta(1)^2)^(1/h0);
    else 
        QUCL=0;
    end     
    S=lamda(X_col-num_pc+1:X_col,X_col-num_pc+1:X_col);
    FAI=P*pinv(S)*P'/T2UCL1+(eye(X_col)-P*P')/QUCL;
    S=cov(Xtrain);
    g=trace((S*FAI)^2)/trace(S*FAI);
    h=(trace(S*FAI))^2/trace((S*FAI)^2);
    %综合控制限
    kesi =g*chi2inv(alpha,h);
    P={P};%转化为元祖
    FAI={FAI};%转化为元祖
    lamda={lamda};%转化为元祖
end
