function [sim]  = eros_measure(A,B,weight)
    % compute mean of each column/variable
     [pc1, ~, ~, ~] = princomp(A);
     [pc2, ~, ~, ~] = princomp(B);
     temp = pc1'*pc2;
     for i=1:size(temp,1)
         temp(i,i) = abs(temp(i,i));
     end
     sim = trace(diag(weight)*temp);
end