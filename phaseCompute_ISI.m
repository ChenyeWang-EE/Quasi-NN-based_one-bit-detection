function y = phaseCompute_ISI(x)
len = length(x);
y = zeros(len-1,1);
summ = 0;
for ii = 1:len-1
    if summ == 0 && x(ii) == -1 && x(ii+1) == -1
        y(ii) = 0;
    elseif summ == 0 && x(ii) == -1 && x(ii+1) == 1
        y(ii) = 1;
    elseif summ == 0 && x(ii) == 1 && x(ii+1) == -1
        y(ii) = 2;
    elseif summ == 0 && x(ii) == 1 && x(ii+1) == 1
        y(ii) = 3;
        %------------------------------------------------%
    elseif summ == 1 && x(ii) == -1 && x(ii+1) == -1
        y(ii) = 4;
    elseif summ == 1 && x(ii) == -1 && x(ii+1) == 1
        y(ii) = 5;
    elseif summ == 1 && x(ii) == 1 && x(ii+1) == -1
        y(ii) = 6;
    elseif summ == 1 && x(ii) == 1 && x(ii+1) == 1
        y(ii) = 7;
        %------------------------------------------------%
    elseif summ == 2 && x(ii) == -1 && x(ii+1) == -1
        y(ii) = 8;
    elseif summ == 2 && x(ii) == -1 && x(ii+1) == 1
        y(ii) = 9;
    elseif summ == 2 && x(ii) == 1 && x(ii+1) == -1
        y(ii) = 10;
    elseif summ == 2 && x(ii) == 1 && x(ii+1) == 1
        y(ii) = 11;
        %------------------------------------------------%
    elseif summ == 3 && x(ii) == -1 && x(ii+1) == -1
        y(ii) = 12;
    elseif summ == 3 && x(ii) == -1 && x(ii+1) == 1
        y(ii) = 13;
    elseif summ == 3 && x(ii) == 1 && x(ii+1) == -1
        y(ii) = 14;
    elseif summ == 3 && x(ii) == 1 && x(ii+1) == 1
        y(ii) = 15;
    end
    summ = mod(summ + x(ii),4);
end
end