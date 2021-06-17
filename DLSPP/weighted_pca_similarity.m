function [cost]=weighted_pca_similarity(s1,s2,ndim)
    %计算PCA相似性因子
    [pc1,~,~] = princomp(s1);
    [pc2,~,~] = princomp(s2);
    L = pc1(:,1:ndim);
    M = pc2(:,1:ndim);
    weight = ones(ndim,1);
    weightsum = 0;
    for i=1:ndim
        weight(i) = 1/i;
        weightsum  = weightsum +  1/i;
    end
    A = M'*L;
    for i=1:size(A,1)
        A(i,i) = abs(A(i,i));
    end
    cost=trace(diag(weight/weightsum)*A); 
end
