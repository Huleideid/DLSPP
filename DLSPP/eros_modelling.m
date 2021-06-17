function [eval,evec]  = eros_modelling(x)
    [~, J, K] = size(x);
    eval = zeros(K,J);
    evec = zeros(J,J,K);
    for i = 1:K
        [pc, ~, latent, ~] = princomp(x(:,:,i));
        eval(i,:) = latent.^2;
        evec(:,:,i) = pc; 
    end  
end