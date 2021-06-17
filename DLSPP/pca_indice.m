% function [cost,pc,latent]  = pcaresid(x,ndim,flag);
function [cost,pc,avg]  = pcaresid2(x,ndim)
    [m,n] = size(x);
    if prod(size(ndim)) > 1
        error('Requires a scalar second input.');
    end
    if ndim >= n
        error('The number of columns in the first input must exceed the value of the 2nd input.');
    end
    x(isnan(x))=0;
    avg = mean(x);
    avgx = avg(ones(m,1),:);
    centerx = (x - avg(ones(m,1),:));
    xcov=cov(centerx);
    [~,~,pc] = svd(xcov,0);
    retain = pc(:,1:ndim)';
    predictx = avgx +centerx*(retain'*retain);
    residuals = x - predictx;
    retain=pc(:,ndim+1:end);
    FAI=retain*retain';
    g=trace((xcov*FAI)^2)/(2*trace(xcov*FAI));
    h=2*(trace(xcov*FAI))^2/trace((xcov*FAI)^2);
    alpha=0.90;
    %×ÛºÏ¿ØÖÆÏÞ
    if isnan(g)==0
        cost =g*chi2inv(alpha,h);
    else
        cost=0;
    end
end
% mean(sum(residuals.^2,2))+mean(sum(score(:,1:ndim).^2,2))=mean(sum( x.^2,2))

