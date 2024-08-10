function y = QAM_label_compute_L3(x,M)
% M : modulation order
%
bitsPerSym = log2(M);
xLen = floor( (length(x)-6) / bitsPerSym +1);
y = zeros(xLen,1);
% mapping = reshape(x, 2*bitsPerSym);
for ii = 1 : xLen
    tmp = x( (ii-1)*bitsPerSym + (1:6 ) ).';
    win1 = [ tmp(3:4),tmp(1:2) ];
    win2 = [ tmp(5:6),tmp(3:4) ];
    win = [win1,win2];
    tmp = replace(num2str(win),' ','');
    res = str2double(tmp);
    switch res
        % ------------------------%
        case 00000000
            y(ii) = 0;
        case 00000100
            y(ii) = 1;
        case 00001000
            y(ii) = 2;
        case 00001100
            y(ii) = 3;
        % ------------------------%
        case 00010000
            y(ii) = 4;
        case 00010100
            y(ii) = 5;
        case 00011000
            y(ii) = 6;
        case 00011100
            y(ii) = 7;
        % ------------------------%
        case 00100000
            y(ii) = 8;
        case 00100100
            y(ii) = 9;
        case 00101000
            y(ii) = 10;
        case 00101100
            y(ii) = 11;

        % ------------------------%
        case 00110000
            y(ii) = 12;
        case 00110100
            y(ii) = 13;
        case 00111000
            y(ii) = 14;
        case 00111100
            y(ii) = 15;
        % ------------------------%
        case 01000001
            y(ii) = 16;
        case 01000101
            y(ii) = 17;
        case 01001001
            y(ii) = 18;
        case 01001101
            y(ii) = 19;
        % ------------------------%
        case 01010001
            y(ii) = 20;
        case 01010101
            y(ii) = 21;
        case 01011001
            y(ii) = 22;
        case 01011101
            y(ii) = 23;
        % ------------------------%
        case 01100001
            y(ii) = 24;
        case 01100101
            y(ii) = 25;
        case 01101001
            y(ii) = 26;
        case 01101101
            y(ii) = 27;

        % ------------------------%
        case 01110001
            y(ii) = 28;
        case 01110101
            y(ii) = 29;
        case 01111001
            y(ii) = 30;
        case 01111101
            y(ii) = 31;
        % ------------------------%
        case 10000010
            y(ii) = 32;
        case 10000110
            y(ii) = 33;
        case 10001010
            y(ii) = 34;
        case 10001110
            y(ii) = 35;
        % ------------------------%
        case 10010010
            y(ii) = 36;
        case 10010110
            y(ii) = 37;
        case 10011010
            y(ii) = 38;
        case 10011110
            y(ii) = 39;
        % ------------------------%
        case 10100010
            y(ii) = 40;
        case 10100110
            y(ii) = 41;
        case 10101010
            y(ii) = 42;
        case 10101110
            y(ii) = 43;

        % ------------------------%
        case 10110010
            y(ii) = 44;
        case 10110110
            y(ii) = 45;
        case 10111010
            y(ii) = 46;
        case 10111110
            y(ii) = 47;

        % ------------------------%
        case 11000011
            y(ii) = 48;
        case 11000111
            y(ii) = 49;
        case 11001011
            y(ii) = 50;
        case 11001111
            y(ii) = 51;
        % ------------------------%
        case 11010011
            y(ii) = 52;
        case 11010111
            y(ii) = 53;
        case 11011011
            y(ii) = 54;
        case 11011111
            y(ii) = 55;
        % ------------------------%
        case 11100011
            y(ii) = 56;
        case 11100111
            y(ii) = 57;
        case 11101011
            y(ii) = 58;
        case 11101111
            y(ii) = 59;

        % ------------------------%
        case 11110011
            y(ii) = 60;
        case 11110111
            y(ii) = 61;
        case 11111011
            y(ii) = 62;
        case 11111111
            y(ii) = 63;

    end
end