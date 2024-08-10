function y = state_tran_output_4QAM_L3(hISI)
    y = zeros(16,4);
    map= [-1 / sqrt(2) + 1i / sqrt(2);
          -1 / sqrt(2) - 1i / sqrt(2);
           1 / sqrt(2) + 1i / sqrt(2);
           1 / sqrt(2) - 1i / sqrt(2)];

    y(1,1) = hISI(1) * map(1) + hISI(2) * map(1) + hISI(3) * map(1);
    y(1,2) = hISI(1) * map(2) + hISI(2) * map(1) + hISI(3) * map(1);
    y(1,3) = hISI(1) * map(3) + hISI(2) * map(1) + hISI(3) * map(1);
    y(1,4) = hISI(1) * map(4) + hISI(2) * map(1) + hISI(3) * map(1);

    y(2,1) = hISI(1) * map(1) + hISI(2) * map(1) + hISI(3) * map(2);
    y(2,2) = hISI(1) * map(2) + hISI(2) * map(1) + hISI(3) * map(2);
    y(2,3) = hISI(1) * map(3) + hISI(2) * map(1) + hISI(3) * map(2);
    y(2,4) = hISI(1) * map(4) + hISI(2) * map(1) + hISI(3) * map(2);

    y(3,1) = hISI(1) * map(1) + hISI(2) * map(1) + hISI(3) * map(3);
    y(3,2) = hISI(1) * map(2) + hISI(2) * map(1) + hISI(3) * map(3);
    y(3,3) = hISI(1) * map(3) + hISI(2) * map(1) + hISI(3) * map(3);
    y(3,4) = hISI(1) * map(4) + hISI(2) * map(1) + hISI(3) * map(3);

    y(4,1) = hISI(1) * map(1) + hISI(2) * map(1) + hISI(3) * map(4);
    y(4,2) = hISI(1) * map(2) + hISI(2) * map(1) + hISI(3) * map(4);
    y(4,3) = hISI(1) * map(3) + hISI(2) * map(1) + hISI(3) * map(4);
    y(4,4) = hISI(1) * map(4) + hISI(2) * map(1) + hISI(3) * map(4);

    y(5,1) = hISI(1) * map(1) + hISI(2) * map(2) + hISI(3) * map(1);
    y(5,2) = hISI(1) * map(2) + hISI(2) * map(2) + hISI(3) * map(1);
    y(5,3) = hISI(1) * map(3) + hISI(2) * map(2) + hISI(3) * map(1);
    y(5,4) = hISI(1) * map(4) + hISI(2) * map(2) + hISI(3) * map(1);

    y(6,1) = hISI(1) * map(1) + hISI(2) * map(2) + hISI(3) * map(2);
    y(6,2) = hISI(1) * map(2) + hISI(2) * map(2) + hISI(3) * map(2);
    y(6,3) = hISI(1) * map(3) + hISI(2) * map(2) + hISI(3) * map(2);
    y(6,4) = hISI(1) * map(4) + hISI(2) * map(2) + hISI(3) * map(2);

    y(7,1) = hISI(1) * map(1) + hISI(2) * map(2) + hISI(3) * map(3);
    y(7,2) = hISI(1) * map(2) + hISI(2) * map(2) + hISI(3) * map(3);
    y(7,3) = hISI(1) * map(3) + hISI(2) * map(2) + hISI(3) * map(3);
    y(7,4) = hISI(1) * map(4) + hISI(2) * map(2) + hISI(3) * map(3);

    y(8,1) = hISI(1) * map(1) + hISI(2) * map(2) + hISI(3) * map(4);
    y(8,2) = hISI(1) * map(2) + hISI(2) * map(2) + hISI(3) * map(4);
    y(8,3) = hISI(1) * map(3) + hISI(2) * map(2) + hISI(3) * map(4);
    y(8,4) = hISI(1) * map(4) + hISI(2) * map(2) + hISI(3) * map(4);

    y(9,1) = hISI(1) * map(1) + hISI(2) * map(3) + hISI(3) * map(1);
    y(9,2) = hISI(1) * map(2) + hISI(2) * map(3) + hISI(3) * map(1);
    y(9,3) = hISI(1) * map(3) + hISI(2) * map(3) + hISI(3) * map(1);
    y(9,4) = hISI(1) * map(4) + hISI(2) * map(3) + hISI(3) * map(1);

    y(10,1) = hISI(1) * map(1) + hISI(2) * map(3) + hISI(3) * map(2);
    y(10,2) = hISI(1) * map(2) + hISI(2) * map(3) + hISI(3) * map(2);
    y(10,3) = hISI(1) * map(3) + hISI(2) * map(3) + hISI(3) * map(2);
    y(10,4) = hISI(1) * map(4) + hISI(2) * map(3) + hISI(3) * map(2);

    y(11,1) = hISI(1) * map(1) + hISI(2) * map(3) + hISI(3) * map(3);
    y(11,2) = hISI(1) * map(2) + hISI(2) * map(3) + hISI(3) * map(3);
    y(11,3) = hISI(1) * map(3) + hISI(2) * map(3) + hISI(3) * map(3);
    y(11,4) = hISI(1) * map(4) + hISI(2) * map(3) + hISI(3) * map(3);

    y(12,1) = hISI(1) * map(1) + hISI(2) * map(3) + hISI(3) * map(4);
    y(12,2) = hISI(1) * map(2) + hISI(2) * map(3) + hISI(3) * map(4);
    y(12,3) = hISI(1) * map(3) + hISI(2) * map(3) + hISI(3) * map(4);
    y(12,4) = hISI(1) * map(4) + hISI(2) * map(3) + hISI(3) * map(4);

    y(13,1) = hISI(1) * map(1) + hISI(2) * map(4) + hISI(3) * map(1);
    y(13,2) = hISI(1) * map(2) + hISI(2) * map(4) + hISI(3) * map(1);
    y(13,3) = hISI(1) * map(3) + hISI(2) * map(4) + hISI(3) * map(1);
    y(13,4) = hISI(1) * map(4) + hISI(2) * map(4) + hISI(3) * map(1);

    y(14,1) = hISI(1) * map(1) + hISI(2) * map(4) + hISI(3) * map(2);
    y(14,2) = hISI(1) * map(2) + hISI(2) * map(4) + hISI(3) * map(2);
    y(14,3) = hISI(1) * map(3) + hISI(2) * map(4) + hISI(3) * map(2);
    y(14,4) = hISI(1) * map(4) + hISI(2) * map(4) + hISI(3) * map(2);

    y(15,1) = hISI(1) * map(1) + hISI(2) * map(4) + hISI(3) * map(3);
    y(15,2) = hISI(1) * map(2) + hISI(2) * map(4) + hISI(3) * map(3);
    y(15,3) = hISI(1) * map(3) + hISI(2) * map(4) + hISI(3) * map(3);
    y(15,4) = hISI(1) * map(4) + hISI(2) * map(4) + hISI(3) * map(3);

    y(16,1) = hISI(1) * map(1) + hISI(2) * map(4) + hISI(3) * map(4);
    y(16,2) = hISI(1) * map(2) + hISI(2) * map(4) + hISI(3) * map(4);
    y(16,3) = hISI(1) * map(3) + hISI(2) * map(4) + hISI(3) * map(4);
    y(16,4) = hISI(1) * map(4) + hISI(2) * map(4) + hISI(3) * map(4);

end
