function [P,FAI,num_pc,lamda,T2UCL1,QUCL,kesi]= pca_batch(Xtrain,threhold,alpha)
    %% �������Ļ��ͱ�׼������Ϊ������Ъ����Ԥ��������ִ���˸ò���
    [X_row,X_col] = size(Xtrain); %��Xtrain�С�����               
    %��Э�������
    sigmaXtrain = cov(Xtrain);
    %��Э���������������ֽ⣬lamdaΪ����ֵ���ɵĶԽ���T����Ϊ��λ��������������lamda�е�����ֵһһ��Ӧ��
    [T,lamda] = eig(sigmaXtrain);                            
                                             
    %ȡ�Խ�Ԫ��(���Ϊһ������)����lamdaֵ�������·�תʹ��Ӵ�С���У���Ԫ������ֵΪ1�����ۼƹ�����С��90%��������Ԫ����
    D = flipud(diag(lamda));                                                                           
    num_pc = 1;                                         
    while sum(D(1:num_pc))/sum(D) < threhold  
        num_pc = num_pc +1;
    end   
    %ȡ��lamda���Ӧ����������
    P = T(:,X_col-num_pc+1:X_col); 
   %% ���������
    %�����Ŷ�Ϊ99%��95%ʱ��T2ͳ�ƿ�����                       
    T2UCL1=num_pc*(X_row-1)*(X_row+1)*finv(alpha,num_pc,X_row - num_pc)/(X_row*(X_row - num_pc));
    %���Ŷ�Ϊ95%��Qͳ�ƿ�����
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
    %�ۺϿ�����
    kesi =g*chi2inv(alpha,h);
    P={P};%ת��ΪԪ��
    FAI={FAI};%ת��ΪԪ��
    lamda={lamda};%ת��ΪԪ��
end
