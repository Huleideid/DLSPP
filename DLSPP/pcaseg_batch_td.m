%--------------------------------修改合并cost下降最快的------------------------------------
function [segment, total_wc] = pcaseg_batch_td(data,num_segments,q)
    % Init 
    segment.lx = 1;
    segment.rx = size(data,3);
    sx=data;
    sx=permute(sx,[2 1 3]);
    sx=reshape(sx,size(data,2),[]);
    segment.c = pcaresid(sx',q,1);
    segment.flag = 0;
    total_wc=[segment.c];
    function [min_mc,best_seg,t_seg1mc,t_seg2mc] = topDownSeg(data,segment,num_pc)
        SS = size(data,3);
        min_mc = Inf;
        t_seg1lx = segment.lx;
        if ( segment.rx - segment.lx)==0
            best_seg=1;
            t_seg1mc = Inf;
            t_seg2mc = Inf;
            return ;
        end
        for iter= segment.lx:( segment.rx - 1)%1 to (SS - 1)
            sx1=data(:,:,t_seg1lx:iter);%1 to (SS - 1)
            sx1=permute(sx1,[2 1 3]);
            sx1=reshape(sx1,size(data,2),[]);
            seg1mc = pcaresid(sx1',num_pc,1);
            t_seg2lx = iter + 1;
            t_seg2rx = segment.rx ;          
            sx2=data(:,:,t_seg2lx :t_seg2rx);
            sx2=permute(sx2,[2 1 3]);
            sx2=reshape(sx2,size(data,2),[]); 
            seg2mc = pcaresid(sx2',num_pc,1);
            t_mc = (iter-t_seg1lx+1)*seg1mc + (t_seg2rx-t_seg2lx+1)*seg2mc - segment.c*(t_seg2rx-t_seg1lx+1);
            if(t_mc < min_mc)
                min_mc = t_mc;
                best_seg = iter;
                t_seg1mc = seg1mc;
                t_seg2mc = seg2mc;
            end
        end
    end  %end of function seg = topDownSeg(data,f_lx)
    while  (length(segment) < num_segments) %分割终止条件
        for i=1:length(segment)
            if [segment(i).flag]==0
                [segment(i).min_mc,segment(i).bestseg,segment(i).t_seg1mc,segment(i).t_seg2mc] = topDownSeg(data,segment(i),q);
                segment(i).flag = 1;%记录子序列是否已经寻找最佳分割点    
            end
            if min([segment(:).flag])==1%判断是否所有的子序列均已经寻找最佳分割点
                index=find([segment(:).min_mc]==min([segment(:).min_mc]));
                segment1.lx=segment(index).lx;
                segment1.rx=segment(index).bestseg;
                segment1.t_seg1mc=inf;
                segment1.t_seg2mc=inf;
                segment1.flag=0;
                segment1.c=segment(index).t_seg1mc;
                segment1.min_mc=inf;
                segment1.bestseg=1;
                segment2.lx=segment(index).bestseg+1;
                segment2.rx=segment(index).rx;
                segment2.t_seg1mc=inf;
                segment2.t_seg2mc=inf;
                segment2.flag=0;
                segment2.min_mc=inf;
                segment2.c=segment(index).t_seg2mc;
                segment2.bestseg=1;
                segment=[segment;segment1;segment2];
                segment(index)=[];
                wc=0;
                for j=1:length(segment)
                    wc=wc+(segment(j).rx-segment(j).lx+1)/size(data,3)*segment(j).c;
                end
                total_wc=[total_wc wc];
            end          
        end
    end
    % 最后重新按照子序列进行排序
    temp=[segment(:).lx];
    [~,I] = sort(temp);
    segment=segment(I);
end