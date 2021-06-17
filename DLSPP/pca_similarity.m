function [cost]=pca_similarity(s1,s2,ndim)
    %计算PCA相似性因子
    [pc1,~,~] = princomp(s1);
    [pc2,~,~] = princomp(s2);
    L = pc1(:,1:ndim);
    M = pc2(:,1:ndim);
    cost=trace(L'*M*M'*L)/ndim; 
end
