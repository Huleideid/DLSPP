function [cost,PC,avg,latent]  = pcaresid(x,ndim,flag);
[m,n] = size(x);
if prod(size(ndim)) > 1
    error('Requires a scalar second input.');
end
if ndim >= n
    error('The number of columns in the first input must exceed the value of the 2nd input.');
end
[pc,score,latent] = princomp(x);
avg = mean(x);
avgx = avg(ones(m,1),:);
PC=pc(:,1:ndim);
retain = pc(:,1:ndim)';
predictx = avgx +score(:,1:ndim)*retain;
residuals = x - predictx;
if flag==1 %Q model error
      cost=mean(sum(residuals.^2,2)); 
else       %Hotelling T2
%      cost=mean(mean([score(:,1:ndim).^2]')');
      cost=mean(sum(score(:,1:ndim).^2,2));
end   
% mean(sum(residuals.^2,2))+mean(sum(score(:,1:ndim).^2,2))=mean(sum( x.^2,2))

