% function   [X_mean,X_std,P,FAI,num_pc,lamda,T2UCL1,QUCL,kesi]=subpca_train(X_train,threhold,alpha,lx,rx)
function   [P,FAI,num_pc,lamda,T2UCL1,QUCL,kesi]=multiphase_pca(X_train,threhold,alpha,lx,rx)
% multiphase pca modelling
    K=length(lx);
    [~,J,~]=size(X_train);
    X_std=zeros(K,J);
    X_mean=zeros(K,J);
    T2UCL1=zeros(K,1);
    QUCL=zeros(K,1);
    num_pc=zeros(K,1);
    kesi=zeros(K,1);
    P=cell(K,1);
    FAI=cell(K,1);
    lamda=cell(K,1);
    for i=1:K
        temp=X_train(:,:,lx(i):rx(i));
        temp=permute(temp,[2 1 3]);
        temp=reshape(temp,size(X_train,2),[]);
        [P(i),FAI(i),num_pc(i),lamda(i),T2UCL1(i),QUCL(i),kesi(i)]= pca_batch(temp',threhold,alpha);
%       [X_mean(i,:),X_std(i,:),P(i),FAI(i),num_pc(i),lamda(i),T2UCL1(i),QUCL(i),kesi(i)]= subpca(temp',threhold,alpha);
    end
end    

