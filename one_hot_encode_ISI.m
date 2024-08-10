function [yi] = one_hot_encode_ISI(index)

yi = zeros(16,1);
yi(index+1) = 1; 
end

