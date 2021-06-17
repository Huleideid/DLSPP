function [segment,tc,wc] = pcaseg_batch_bu(data,num_segments,q,flag,min_length)
    % The basic algorithm works by creating
    % a fine segmented representation then merging 
    % the lowest cost segments until only 'num_segments' remain.

    % [segment,tc] = pcaseg(data,num_segments,q,flag)
    % pcaseg.m - Bottom-up segmentation of multivariate time series using PCA based cost functions
    % Created by Zolt¨¢n Bank¨®, 2006
    % Inputs:
    %         data: three-dimension ,the third deminseion is time
    %         num_segments: desired number of segments
    %         q: number of principal compnents to be retained
    %         flag: cost type (0: Hotelling (T2), 1: avarage residual error(Q))
    % Outputs:
    %         segment: a szegment¨¢l¨¢st le¨ªr¨® strukture, elemei a szegmensek
    %         tc: total cost of the resulted segmentation
    %         wc: weighted cost of the resulted segmentation
    if nargin==4
        min_length=1;
    end
    minres=floor(size(data,3)/size(data,3)); %minimal resolution;initial segmentation -> MUST BE MODIFIED BASED ON YOUR DATA!!!
    left_x       = [1 : minres : size(data,3)];    % Find the left x values vector for the "fine segmented representation".
    right_x      = left_x + minres-1;                  % Find the right x values vector for the "fine segmented representation".
    right_x(end) = size(data,3);                     % Special case, the rightmost endpoint.
    number_of_segments = length(left_x );            % Calculate the number of segments in the initial "fine segmented representation".
    if flag==1  % Q based related cost
        % Initialize the segments in the "fine segmented representation". 
        for i = 1 : number_of_segments 
           segment(i).lx = left_x(i);
           segment(i).rx = right_x(i);
           segment(i).mc = inf;
           segment(i).c = inf;
        end;
        tc=[];
        wc=[];

        % Initialize the merge cost of the segments in the "fine segmented representation". 
        % compute merge costs (i.e. cost of two consecutive segments)
        for i = 1 : number_of_segments
           sx=data(:,:,segment(i).lx :segment(i).rx);
           sx=permute(sx,[2 1 3]);
           sx=reshape(sx,size(data,2),[]);
           segment(i).c = pcaresid(sx',q,flag);  % cost of the segment itself
        end
        for i = 1 : number_of_segments -1
            sx=data(:,:,segment(i).lx :segment(i+1).rx);
            sx=permute(sx,[2 1 3]);
            sx=reshape(sx,size(data,2),[]);
            segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
        end
        % Keep merging the lowest cost segments until only 'num_segments' remain. 
        while  length(segment) > num_segments  
           [~, i ] = min([segment(:).mc]);                              % Find the location "i", of the cheapest merge.
           if i > 1 & i < length(segment) -1								% The typical case, neither of the two segments to be merged are end segments
                   segment(i).rx = segment(i+1).rx;
                   segment(i+1) = [];
                   sx=data(:,:,segment(i).lx :segment(i).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).c = pcaresid(sx',q,flag);
                   sx=data(:,:,segment(i).lx :segment(i+1).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx +1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c;                   
                   i = i - 1;
                   sx=data(:,:,segment(i).lx :segment(i+1).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
           elseif i == 1    % Special case: The leftmost segment must be merged.
                   segment(i).rx = segment(i+1).rx;
                   segment(i+1) = [];
                   sx=data(:,:,segment(i).lx :segment(i).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).c = pcaresid(sx',q,flag);
                   sx=data(:,:,segment(i).lx :segment(i+1).rx);
                   sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                   segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
       
           else % special case 3: last segment is merged                                                             % Special case: The rightmost segment must be merged.
                  segment(i).rx = segment(i+1).rx;% update the last data point of the newly merged segment
                  segment(i+1) = []; % delete the second segments of two we just merged
                  sx=data(:,:,segment(i).lx :segment(i).rx);
                  sx=permute(sx,[2 1 3]);
                  sx=reshape(sx,size(data,2),[]);
                  segment(i).c = pcaresid(sx',q,flag);
                  segment(i).mc = inf;
                  i = i - 1;% decrease index
                  sx=data(:,:,segment(i).lx :segment(i+1).rx);
                  sx=permute(sx,[2 1 3]);
                  sx=reshape(sx,size(data,2),[]);
                  segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
           end; % end of if i > 1 & i < length(segment) -1
           % total cost 
           tc=[tc; sum([segment.c])];
           % weighted cost
           segcost=[];
           norm_seglength=[];
           seglength=[];
           for j=1:length(segment) 
                seglength(j) = segment(j).rx - segment(j).lx+1;
                norm_seglength(j) = seglength(j)/right_x(end);
                segcost(j) = segment(j).c;
           end
           wc = [wc; sum([norm_seglength .* segcost])];
        end; % end of while  length(segment) > num_segments
     %% postprepossing step
      while min(seglength)<min_length
            if length(segment)==2
                segment(0).rx = segment(1).rx ;
                segment(1) = [];
                sx=permute(data,[2 1 3]);
                sx=reshape(sx,size(data,2),[]);
                segment(0).c = pcaresid(sx',q,flag); 
                segment(0).mc = inf;
            else % segment must be greatter than 2
                ind=find(seglength < min_length);
                i=ind(1);
                if (i >1) | (i<length(segment))
                    [~, shift] = min([segment(i-1:i).mc]);
                    i = i-2+shift;
                end       
                if  i > 1 & i < length(segment) -1	% The typical case, neither of the two segments to be merged are end segments
                    segment(i).rx = segment(i+1).rx;
                    segment(i+1) = [];
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag);
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag); 
                    i = i - 1;
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag); 
                elseif i==1   % case : first segment 
                    segment(i).rx = segment(i+1).rx;
                    segment(i+1) = [];
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag); 
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag); 
                elseif (i==length(segment) -1) % case : last segment-1                                                          % Special case: The rightmost segment must be merged.
                    segment(i).rx = segment(i+1).rx;% update the last data point of the newly merged segment
                    segment(i+1) = []; % delete the second segments of two we just merged
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag);
                    segment(i).mc = inf;
                    i = i - 1;% decrease index
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
                elseif i==length(segment) %case : last segment 
                    i = i-1;
                    segment(i).rx = segment(i+1).rx;% update the last data point of the newly merged segment
                    segment(i+1) = []; % delete the second segments of two we just merged
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag);
                    segment(i).mc=inf;
                    i=i-1;
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
                end
            end
            % update seglength
            seglength=[];
            for j=1:length(segment) 
                seglength(j) = segment(j).rx - segment(j).lx+1;
            end
       end
    else % T2 based related cost
        % Initialize the segments in the "fine segmented representation". 
        for i = 1 : number_of_segments 
           segment(i).lx = left_x(i);
           segment(i).rx = right_x(i);
           segment(i).mc = -inf;
           segment(i).c = -inf;
        end;
        tc=[];
        wc=[];
        % Initialize the merge cost of the segments in the "fine segmented representation". 
        for i = 1 : number_of_segments
           sx=data(:,:,segment(i).lx :segment(i).rx);
           sx=permute(sx,[2 1 3]);
           sx=reshape(sx,size(data,2),[]);
           segment(i).c = pcaresid(sx',q,flag);  % cost of the segment itself
        end
        % compute merge costs (i.e. cost of two consecutive segments)
        for i = 1 : number_of_segments -1
            sx=data(:,:,segment(i).lx :segment(i+1).rx);
            sx=permute(sx,[2 1 3]);
            sx=reshape(sx,size(data,2),[]);
            segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
        end
        % Keep merging the lowest cost segments until only 'num_segments' remain. 
        while  length(segment) > num_segments  
           [~, i] = max([segment(:).mc]);                              % Find the location "i", of the cheapest merge.
           if i > 1 & i < length(segment) -1								% The typical case, neither of the two segments to be merged are end segments
                   segment(i).rx = segment(i+1).rx;
                   segment(i+1) = [];
                   sx=data(:,:,segment(i).lx :segment(i).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).c = pcaresid(sx',q,flag);
                   sx=data(:,:,segment(i).lx :segment(i+1).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx +1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i+1).rx-segment(i).lx+1)*segment(i).c; 
                   i = i - 1;
                   sx=data(:,:,segment(i).lx :segment(i+1).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
           elseif i == 1    % Special case: The leftmost segment must be merged.
                   segment(i).rx = segment(i+1).rx;
                   segment(i+1) = [];
                   sx=data(:,:,segment(i).lx :segment(i).rx);
                   sx=permute(sx,[2 1 3]);
                   sx=reshape(sx,size(data,2),[]);
                   segment(i).c = pcaresid(sx',q,flag);
                   sx=data(:,:,segment(i).lx :segment(i+1).rx);
                   sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                   segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
           else % special case 3: last segment is merged                                                             % Special case: The rightmost segment must be merged.
                  segment(i).rx = segment(i+1).rx;% update the last data point of the newly merged segment
                  segment(i+1) = []; % delete the second segments of two we just merged
                  sx=data(:,:,segment(i).lx :segment(i).rx);
                  sx=permute(sx,[2 1 3]);
                  sx=reshape(sx,size(data,2),[]);
                  segment(i).c = pcaresid(sx',q,flag);
                  segment(i).mc = -inf;
                  i = i - 1;% decrease index
                  sx=data(:,:,segment(i).lx :segment(i+1).rx);
                  sx=permute(sx,[2 1 3]);
                  sx=reshape(sx,size(data,2),[]);
                  segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
           end; % end of if i > 1 & i < length(segment) -1
           % total cost 
           tc=[tc; sum([segment.c])];
           % weighted cost
           segcost=[];
           norm_seglength=[];
           seglength=[];
           for j=1:length(segment) 
                seglength(j) = segment(j).rx - segment(j).lx+1;
                norm_seglength(j) = seglength(j)/right_x(end);
                segcost(j) = segment(j).c;
           end
           wc = [wc; sum([norm_seglength .* segcost])];
       end; % end of while  length(segment) > num_segments
       while min(seglength)<min_length
            if length(segment)==2
                segment(0).rx = segment(1).rx ;
                segment(1) = [];
                sx=permute(data,[2 1 3]);
                sx=reshape(sx,size(data,2),[]);
                segment(0).c = pcaresid(sx',q,flag); 
                segment(0).mc = -inf;
            else % segment must be greatter than 2
                ind=find(seglength < min_length);
                i=ind(1);
                if (i >1) | (i<length(segment))
                    [~, shift] = max([segment(i-1:i).mc]);
                    i = i-2+shift;
                end     
                if  i > 1 & i < length(segment) -1	% The typical case, neither of the two segments to be merged are end segments
                    segment(i).rx = segment(i+1).rx;
                    segment(i+1) = [];
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag);
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag); 
                    i = i - 1;
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag); 
                elseif i==1   % case : first segment 
                    segment(i).rx = segment(i+1).rx;
                    segment(i+1) = [];
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag); 
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag); 
                elseif (i==length(segment) -1) % case : last  last segment                                                        % Special case: The rightmost segment must be merged.
                    segment(i).rx = segment(i+1).rx;% update the last data point of the newly merged segment
                    segment(i+1) = []; % delete the second segments of two we just merged
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag);
                    segment(i).mc = -inf;
                    i = i - 1;% decrease index
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
                elseif i==length(segment) %case : last segment 
                    i = i-1;
                    segment(i).rx = segment(i+1).rx;% update the last data point of the newly merged segment
                    segment(i+1) = []; % delete the second segments of two we just merged
                    sx=data(:,:,segment(i).lx :segment(i).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).c = pcaresid(sx',q,flag);
                    segment(i).mc=-inf;
                    i=i-1;
                    sx=data(:,:,segment(i).lx :segment(i+1).rx);
                    sx=permute(sx,[2 1 3]);
                    sx=reshape(sx,size(data,2),[]);
                    segment(i).mc = pcaresid(sx',q,flag)*(segment(i+1).rx-segment(i).lx+1)-(segment(i+1).rx-segment(i+1).lx+1)*segment(i+1).c-(segment(i).rx-segment(i).lx+1)*segment(i).c; % compute merge cost with the consecutive segment
                end
            end
            % update seglength
            seglength=[];
            for j=1:length(segment) 
                seglength(j) = segment(j).rx - segment(j).lx+1;
            end
       end
    end% end of if flag==1
end

