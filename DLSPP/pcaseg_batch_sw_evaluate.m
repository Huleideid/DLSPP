function [wc] = pcaseg_batch_sw_evaluate(test_data,q,train_segment)
    %First init 
    flag=1;
    for i=1:length(train_segment)
        sx=test_data(:,:,train_segment(i).lx :train_segment(i).rx);
        sx=permute(sx,[2 1 3]);
        sx=reshape(sx,size(test_data,2),[]); 
        test_segment(i).c = pcaresid(sx',q,flag);  % cost of the train_segment itself
    end
   % Global cost
   segcost=[];
   norm_seglength=[];
   seglength=[];
   for j=1:length(test_segment) 
        seglength(j) = train_segment(j).rx - train_segment(j).lx+1;
        norm_seglength(j) = seglength(j)/(size(test_data,3));
        segcost(j) = test_segment(j).c;
   end  
   wc=sum([(norm_seglength).* segcost]);
end