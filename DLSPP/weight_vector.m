
%% compute weight vector from eigen values
function [normed_weight_vector]  = weight_vector(eval)
    % compute mean of each column/variable
    weight = mean(eval,1);
    % normalize weight vector
    total_weight = sum(weight);
    normed_weight_vector = weight/total_weight;
end