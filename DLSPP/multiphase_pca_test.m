% function   [T2,Q]=subpca_test(X_test,X_mean,X_std,P,num_pc,lamda,lx,rx)% 训练
% function   [T2,Q]=subpca_test(X_test,P,num_pc,lamda,lx,rx)% 训练
function   [T2,Q,fai]=multiphase_pca_test(X_test,P,num_pc, lamda, T2UCL, QUCL, lx, rx)% 训练
    [num,J,K]=size(X_test);
    T2=zeros(num,K);
    Q=zeros(num,K);
    fai=zeros(num,K);
    c=length(lx);
    for i=1:c
    %求T2统计量，Q统计量
        for j=lx(i):rx(i)
            for m=1:num
                [r,y] = size(P{i}*P{i}');
                I = eye(r,y);
                T2(m,j)=X_test(m,:,j)* P{i}*pinv(lamda{i}((J-num_pc(i)+1):J,J-num_pc(i)+1:J))* P{i}'*X_test(m,:,j)';  
                Q(m,j) = X_test(m,:,j)*(I - P{i}*P{i}')*(I - P{i}*P{i}')'*X_test(m,:,j)';
                fai(m,j)=(Q(m,j)/QUCL(1,j))+(T2(m,j)/T2UCL(1,j));
            end
        end
    end
end


