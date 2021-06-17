function [stdx,avg,var] = zscore_bactch(x,method) 
%-------------------------------------------------------------------------- 
    switch method
        case 1 %方式2
            temp=permute(x,[2,3,1]);
            p=reshape(temp,size(temp,1),[]);
            temp1=zscore(p');
            temp1=temp1';
            finalp=reshape(temp1,size(temp,1),size(temp,2),size(temp,3));
            finalp=permute(finalp,[1,3,2]);
            stdx=permute(finalp,[2,1,3]);
        case 2  %方式2
            for i=1:size(x,3) 
                avg(i,:)=mean(x(:,:,i));          
                var(i,:)=std(x(:,:,i));
                for j=1:size(x,2)   
                % if std(x(:,j,i))<=1*10^(-6)
                    if std(x(:,j,i))<=1*10^(-10)
                        stdx(:,j,i)=x(:,j,i)-mean(x(:,j,i));
                    else
                        stdx(:,j,i)=zscore(x(:,j,i));
                    end
                end
            end
            if nargout < 2, return; end
            for i=1:size(x,3) 
                avg(i,:)=mean(x(:,:,i));
                var(i,:)=std(x(:,:,i));
             end
    end
end
