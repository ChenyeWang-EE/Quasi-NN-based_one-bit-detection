function [yi] = one_hot_encode_ISI_L3(index)

yi = zeros(64,1);
yi(index+1) = 1; 
end

