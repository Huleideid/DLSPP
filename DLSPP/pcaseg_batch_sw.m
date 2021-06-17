function [segment, wc] = pcaseg_batch_sw(data,minlength,q,alpha)
    % �������ڷ���ʹ��speͳ����
    %% Init 
    t_lx = 1;
    segment.lx=t_lx;
    segment.rx=0;
    segment.mc=inf;
    function seg = slideWindowSeg(data,num_pc,f_lx)
        SS = size(data,3);
        iter = minlength;%2 to SS - 1
        sx=data(:,:,1:iter);
        sx=permute(sx,[2 1 3]);
        sx=reshape(sx,size(data,2),[]);
        t_mc = pca_indice(sx',num_pc);
        count=0;
        for iter=(minlength+1):(SS - minlength)
            seg.lx = 1;
            seg.rx = iter;
            sx1=data(:,:,seg.lx :seg.rx);
            sx1=permute(sx1,[2 1 3]);
            sx1=reshape(sx1,size(data,2),[]);
            seg.mc = pca_indice(sx1',num_pc);
            if(seg.mc  >= alpha*t_mc)
                count=count+1;
                if count==3
                    seg.rx =seg.rx -3;
                    break;
                end 
            else
                count=0; %reset flag
            end
        end  %end of iter=(minlength+1):(SS - minlength)
        if iter==(SS - minlength)
            seg.rx=SS;
        end
        seg.lx = f_lx;%����ƫ��
        seg.rx = seg.rx + f_lx - 1;%����ƫ��
    end  %end of function seg = topDownSeg(data,f_lx)
    while  (size(data,3)-segment(end).rx) >= (2*minlength+3)%������������������Ŷ����н��зָ�
        new_seg = slideWindowSeg(data(:,:,t_lx :end),q,t_lx);
        tc=new_seg.mc;
        segment = [segment new_seg];
        t_lx  = new_seg.rx+1; 
    end
    if t_lx~=(size(data,3)+1)  %����������и���,ֱ��ȡ�յ����������
        sx=data(:,:,t_lx:end);
        sx=permute(sx,[2 1 3]);
        sx=reshape(sx,size(data,2),[]);
        new_seg.mc=pca_indice(sx',q);
        new_seg.lx=t_lx;
        new_seg.rx=size(data,3);
        segment=[segment new_seg];
    end
    segment(1)=[];
    %% ����ȫ�ִ��ۺ���  calculate Global cost
    for i=1:length(segment)
        sx=data(:,:,segment(i).lx :segment(i).rx);
        sx=permute(sx,[2 1 3]);
        sx=reshape(sx,size(data,2),[]); 
        segment(i).c = pcaresid(sx',q,1);  % cost of the train_segment itself
    end
    segcost=[];
    norm_seglength=[];
    seglength=[];
    for j=1:length(segment) 
        seglength(j) = segment(j).rx - segment(j).lx+1;
        norm_seglength(j) = seglength(j)/(size(data,3));
        segcost(j) = segment(j).c;
    end  
    wc=sum([(norm_seglength).* segcost]);
end