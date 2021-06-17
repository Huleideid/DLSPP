function [test_segment,wc] = pcaseg_batch_td_evaluate(data,num_segments,q,train_segment,test_data)
    %First init 
    wc=[];
    flag=1;
    for i=1:length(train_segment)
        sx=test_data(:,:,train_segment(i).lx :train_segment(i).rx);
        sx=permute(sx,[2 1 3]);
        sx=reshape(sx,size(test_data,2),[]); 
        test_segment(i).c = pcaresid(sx',q,flag);  % cost of the train_segment itself
    end
   % weighted cost
   segcost=[];
   norm_seglength=[];
   seglength=[];
   for j=1:length(test_segment) 
        seglength(j) = train_segment(j).rx - train_segment(j).lx+1;
        norm_seglength(j) = seglength(j)/(size(test_data,3));
        segcost(j) = test_segment(j).c;
   end  
   wc = [wc; sum([(norm_seglength).* segcost])];
    function [min_mc,best_seg,t_seg1mc,t_seg2mc] = topDownSeg(data,train_segment,num_pc)
        SS = size(data,3);
        min_mc = Inf;
        t_seg1lx = train_segment.lx;
        if ( train_segment.rx - train_segment.lx)==0
            best_seg=1;
            t_seg1mc = Inf;
            t_seg2mc = Inf;
            return ;
        end
        for iter= train_segment.lx:( train_segment.rx - 1)%2 to SS - 1
            sx1=data(:,:,t_seg1lx:iter);%2 to SS - 1)
            sx1=permute(sx1,[2 1 3]);
            sx1=reshape(sx1,size(data,2),[]);
            seg1mc = pcaresid(sx1',num_pc,1);
            t_seg2lx = iter + 1;
            t_seg2rx = train_segment.rx ;          
            sx2=data(:,:,t_seg2lx :t_seg2rx);
            sx2=permute(sx2,[2 1 3]);
            sx2=reshape(sx2,size(data,2),[]); 
            seg2mc = pcaresid(sx2',num_pc,1);
            t_mc = (iter-t_seg1lx+1)*seg1mc + (t_seg2rx-t_seg2lx+1)*seg2mc-train_segment.c*(t_seg2rx-t_seg1lx+1);
            if(t_mc<min_mc)
                min_mc =t_mc;
                best_seg=iter;
                t_seg1mc=seg1mc;
                t_seg2mc=seg2mc;
            end
        end
    end  %end of function seg = topDownSeg(data,f_lx)
    while  (length(train_segment) < num_segments)
        for i=1:length(train_segment)
            if train_segment(i).flag==0
                [train_segment(i).min_mc,train_segment(i).bestseg,train_segment(i).t_seg1mc,train_segment(i).t_seg2mc] = topDownSeg(data,train_segment(i),q);
                train_segment(i).flag = 1;
            end
            if min([train_segment(:).flag])==1
                index=find([train_segment(:).min_mc]==min([train_segment(:).min_mc]));
                segment1.lx=train_segment(index).lx;
                segment1.rx=train_segment(index).bestseg;
                segment1.t_seg1mc=inf;
                segment1.t_seg2mc=inf;
                segment1.flag=0;
                segment1.c=train_segment(index).t_seg1mc;
                segment1.min_mc=inf;
                segment1.bestseg=1;
                segment2.lx=train_segment(index).bestseg+1;
                segment2.rx=train_segment(index).rx;
                segment2.t_seg1mc=inf;
                segment2.t_seg2mc=inf;
                segment2.flag=0;
                segment2.min_mc=inf;
                segment2.c=train_segment(index).t_seg2mc;
                segment2.bestseg=1;
                train_segment=[train_segment;segment1;segment2];
                train_segment(index)=[];
                for i=1:length(train_segment)
                    sx=test_data(:,:,train_segment(i).lx :train_segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(test_data,2),[]); 
                    test_segment(i).c = pcaresid(sx',q,flag);  % cost of the train_segment itself
                end
                c=0;
                for j=1:length(train_segment)
                    c=c+(train_segment(j).rx-train_segment(j).lx+1)/size(data,3)*test_segment(j).c;
                end
                wc=[wc c];
            end          
        end
    end
    temp=[train_segment(:).lx];
    [B,I] = sort(temp);
    test_segment=test_segment(I);
end