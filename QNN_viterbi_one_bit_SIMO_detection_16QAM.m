function [DecodeBit,FinalPath] = QNN_viterbi_one_bit_SIMO_detection_16QAM(zkVec,Re_hMtr,Im_hMtr)
%% upsample = 1;
% len = length(zk);
% yprob = zeros(8,8,len);

Path1 = [];
Path2 = [];
Path3 = [];
Path4 = [];
Path5 = [];
Path6 = [];
Path7 = [];
Path8 = [];
Path9 = [];
Path10 = [];
Path11 = [];
Path12 = [];
Path13 = [];
Path14 = [];
Path15 = [];
Path16 = [];

PathMetric = zeros(1,16);
% PathMetric(1) = 0;
% PathMetric(2) = 0;
% numBits = length(zkVec);
numBits = size(zkVec,1);
% % test
BranchMetricMtr = QNN_distance_branch_one_bit_SIMO_16QAM( zkVec,Re_hMtr,Im_hMtr);
for i = 1 : numBits
    NewPathMetric = zeros(16,1);
    % BranchMetric  = distance_branch_one_bit_L2( zk(i),channelH,noiseVar); % distance -log(p)
    % BranchMetric = distance_branch_4QAM( zk,channelH,noiseVar);
    BranchMetric = reshape( BranchMetricMtr(i,:), 16, 16).'; 
    %% viterbi decode
    % update P1
    cost1 = [PathMetric(1) + BranchMetric(1,1);
             PathMetric(2) + BranchMetric(2,1);
             PathMetric(3) + BranchMetric(3,1);
             PathMetric(4) + BranchMetric(4,1);
             PathMetric(5) + BranchMetric(5,1);
             PathMetric(6) + BranchMetric(6,1);
             PathMetric(7) + BranchMetric(7,1);
             PathMetric(8) + BranchMetric(8,1);
             PathMetric(9) + BranchMetric(9,1);
             PathMetric(10) + BranchMetric(10,1);
             PathMetric(11) + BranchMetric(11,1);
             PathMetric(12) + BranchMetric(12,1);
             PathMetric(13) + BranchMetric(13,1);
             PathMetric(14) + BranchMetric(14,1);
             PathMetric(15) + BranchMetric(15,1);
             PathMetric(16) + BranchMetric(16,1)];
    [NewPathMetric(1),idx1] = min( cost1 );
    switch idx1
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP1 = [Path,idx1];
    
    % update P2
    cost2 = [PathMetric(1) + BranchMetric(1,2);
             PathMetric(2) + BranchMetric(2,2);
             PathMetric(3) + BranchMetric(3,2);
             PathMetric(4) + BranchMetric(4,2);
             PathMetric(5) + BranchMetric(5,2);
             PathMetric(6) + BranchMetric(6,2);
             PathMetric(7) + BranchMetric(7,2);
             PathMetric(8) + BranchMetric(8,2);
             PathMetric(9) + BranchMetric(9,2);
             PathMetric(10) + BranchMetric(10,2);
             PathMetric(11) + BranchMetric(11,2);
             PathMetric(12) + BranchMetric(12,2);
             PathMetric(13) + BranchMetric(13,2);
             PathMetric(14) + BranchMetric(14,2);
             PathMetric(15) + BranchMetric(15,2);
             PathMetric(16) + BranchMetric(16,2)];
    [NewPathMetric(2),idx2] = min( cost2 );
    switch idx2
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP2 = [Path,idx2];

    % update P3
    cost3 = [PathMetric(1) + BranchMetric(1,3);
             PathMetric(2) + BranchMetric(2,3);
             PathMetric(3) + BranchMetric(3,3);
             PathMetric(4) + BranchMetric(4,3);
             PathMetric(5) + BranchMetric(5,3);
             PathMetric(6) + BranchMetric(6,3);
             PathMetric(7) + BranchMetric(7,3);
             PathMetric(8) + BranchMetric(8,3);
             PathMetric(9) + BranchMetric(9,3);
             PathMetric(10) + BranchMetric(10,3);
             PathMetric(11) + BranchMetric(11,3);
             PathMetric(12) + BranchMetric(12,3);
             PathMetric(13) + BranchMetric(13,3);
             PathMetric(14) + BranchMetric(14,3);
             PathMetric(15) + BranchMetric(15,3);
             PathMetric(16) + BranchMetric(16,3)];
    [NewPathMetric(3),idx3] = min( cost3 );
    switch idx3
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP3 = [Path,idx3];

    % update P4
    cost4 = [PathMetric(1) + BranchMetric(1,4);
             PathMetric(2) + BranchMetric(2,4);
             PathMetric(3) + BranchMetric(3,4);
             PathMetric(4) + BranchMetric(4,4);
             PathMetric(5) + BranchMetric(5,4);
             PathMetric(6) + BranchMetric(6,4);
             PathMetric(7) + BranchMetric(7,4);
             PathMetric(8) + BranchMetric(8,4);
             PathMetric(9) + BranchMetric(9,4);
             PathMetric(10) + BranchMetric(10,4);
             PathMetric(11) + BranchMetric(11,4);
             PathMetric(12) + BranchMetric(12,4);
             PathMetric(13) + BranchMetric(13,4);
             PathMetric(14) + BranchMetric(14,4);
             PathMetric(15) + BranchMetric(15,4);
             PathMetric(16) + BranchMetric(16,4)];
    [NewPathMetric(4),idx4] = min( cost4 );
    switch idx4
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP4 = [Path,idx4];

    % update P5
    cost5 = [PathMetric(1) + BranchMetric(1,5);
             PathMetric(2) + BranchMetric(2,5);
             PathMetric(3) + BranchMetric(3,5);
             PathMetric(4) + BranchMetric(4,5);
             PathMetric(5) + BranchMetric(5,5);
             PathMetric(6) + BranchMetric(6,5);
             PathMetric(7) + BranchMetric(7,5);
             PathMetric(8) + BranchMetric(8,5);
             PathMetric(9) + BranchMetric(9,5);
             PathMetric(10) + BranchMetric(10,5);
             PathMetric(11) + BranchMetric(11,5);
             PathMetric(12) + BranchMetric(12,5);
             PathMetric(13) + BranchMetric(13,5);
             PathMetric(14) + BranchMetric(14,5);
             PathMetric(15) + BranchMetric(15,5);
             PathMetric(16) + BranchMetric(16,5)];
    [NewPathMetric(5),idx5] = min( cost5 );
    switch idx5
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP5 = [Path,idx5];

    % update P6
    cost6 = [PathMetric(1) + BranchMetric(1,6);
             PathMetric(2) + BranchMetric(2,6);
             PathMetric(3) + BranchMetric(3,6);
             PathMetric(4) + BranchMetric(4,6);
             PathMetric(5) + BranchMetric(5,6);
             PathMetric(6) + BranchMetric(6,6);
             PathMetric(7) + BranchMetric(7,6);
             PathMetric(8) + BranchMetric(8,6);
             PathMetric(9) + BranchMetric(9,6);
             PathMetric(10) + BranchMetric(10,6);
             PathMetric(11) + BranchMetric(11,6);
             PathMetric(12) + BranchMetric(12,6);
             PathMetric(13) + BranchMetric(13,6);
             PathMetric(14) + BranchMetric(14,6);
             PathMetric(15) + BranchMetric(15,6);
             PathMetric(16) + BranchMetric(16,6)];
    [NewPathMetric(6),idx6] = min( cost6 );
    switch idx6
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP6 = [Path,idx6];

    % update P7
    cost7 = [PathMetric(1) + BranchMetric(1,7);
             PathMetric(2) + BranchMetric(2,7);
             PathMetric(3) + BranchMetric(3,7);
             PathMetric(4) + BranchMetric(4,7);
             PathMetric(5) + BranchMetric(5,7);
             PathMetric(6) + BranchMetric(6,7);
             PathMetric(7) + BranchMetric(7,7);
             PathMetric(8) + BranchMetric(8,7);
             PathMetric(9) + BranchMetric(9,7);
             PathMetric(10) + BranchMetric(10,7);
             PathMetric(11) + BranchMetric(11,7);
             PathMetric(12) + BranchMetric(12,7);
             PathMetric(13) + BranchMetric(13,7);
             PathMetric(14) + BranchMetric(14,7);
             PathMetric(15) + BranchMetric(15,7);
             PathMetric(16) + BranchMetric(16,7)];
    [NewPathMetric(7),idx7] = min( cost7 );
    switch idx7
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP7 = [Path,idx7];

    % update P8
    cost8 = [PathMetric(1) + BranchMetric(1,8);
             PathMetric(2) + BranchMetric(2,8);
             PathMetric(3) + BranchMetric(3,8);
             PathMetric(4) + BranchMetric(4,8);
             PathMetric(5) + BranchMetric(5,8);
             PathMetric(6) + BranchMetric(6,8);
             PathMetric(7) + BranchMetric(7,8);
             PathMetric(8) + BranchMetric(8,8);
             PathMetric(9) + BranchMetric(9,8);
             PathMetric(10) + BranchMetric(10,8);
             PathMetric(11) + BranchMetric(11,8);
             PathMetric(12) + BranchMetric(12,8);
             PathMetric(13) + BranchMetric(13,8);
             PathMetric(14) + BranchMetric(14,8);
             PathMetric(15) + BranchMetric(15,8);
             PathMetric(16) + BranchMetric(16,8)];
    [NewPathMetric(8),idx8] = min( cost8 );
    switch idx8
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP8 = [Path,idx8];

    % update P9
    cost9 = [PathMetric(1) + BranchMetric(1,9);
             PathMetric(2) + BranchMetric(2,9);
             PathMetric(3) + BranchMetric(3,9);
             PathMetric(4) + BranchMetric(4,9);
             PathMetric(5) + BranchMetric(5,9);
             PathMetric(6) + BranchMetric(6,9);
             PathMetric(7) + BranchMetric(7,9);
             PathMetric(8) + BranchMetric(8,9);
             PathMetric(9) + BranchMetric(9,9);
             PathMetric(10) + BranchMetric(10,9);
             PathMetric(11) + BranchMetric(11,9);
             PathMetric(12) + BranchMetric(12,9);
             PathMetric(13) + BranchMetric(13,9);
             PathMetric(14) + BranchMetric(14,9);
             PathMetric(15) + BranchMetric(15,9);
             PathMetric(16) + BranchMetric(16,9)];
    [NewPathMetric(9),idx9] = min( cost9 );
    switch idx9
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP9 = [Path,idx9];

    % update P10
    cost10 = [PathMetric(1) + BranchMetric(1,10);
             PathMetric(2) + BranchMetric(2,10);
             PathMetric(3) + BranchMetric(3,10);
             PathMetric(4) + BranchMetric(4,10);
             PathMetric(5) + BranchMetric(5,10);
             PathMetric(6) + BranchMetric(6,10);
             PathMetric(7) + BranchMetric(7,10);
             PathMetric(8) + BranchMetric(8,10);
             PathMetric(9) + BranchMetric(9,10);
             PathMetric(10) + BranchMetric(10,10);
             PathMetric(11) + BranchMetric(11,10);
             PathMetric(12) + BranchMetric(12,10);
             PathMetric(13) + BranchMetric(13,10);
             PathMetric(14) + BranchMetric(14,10);
             PathMetric(15) + BranchMetric(15,10);
             PathMetric(16) + BranchMetric(16,10)];
    [NewPathMetric(10),idx10] = min( cost10 );
    switch idx10
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP10 = [Path,idx10];

    % update P11
    cost11 = [PathMetric(1) + BranchMetric(1,11);
             PathMetric(2) + BranchMetric(2,11);
             PathMetric(3) + BranchMetric(3,11);
             PathMetric(4) + BranchMetric(4,11);
             PathMetric(5) + BranchMetric(5,11);
             PathMetric(6) + BranchMetric(6,11);
             PathMetric(7) + BranchMetric(7,11);
             PathMetric(8) + BranchMetric(8,11);
             PathMetric(9) + BranchMetric(9,11);
             PathMetric(10) + BranchMetric(10,11);
             PathMetric(11) + BranchMetric(11,11);
             PathMetric(12) + BranchMetric(12,11);
             PathMetric(13) + BranchMetric(13,11);
             PathMetric(14) + BranchMetric(14,11);
             PathMetric(15) + BranchMetric(15,11);
             PathMetric(16) + BranchMetric(16,11)];
    [NewPathMetric(11),idx11] = min( cost11 );
    switch idx11
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP11 = [Path,idx11];

    % update P12
    cost12 = [PathMetric(1) + BranchMetric(1,12);
             PathMetric(2) + BranchMetric(2,12);
             PathMetric(3) + BranchMetric(3,12);
             PathMetric(4) + BranchMetric(4,12);
             PathMetric(5) + BranchMetric(5,12);
             PathMetric(6) + BranchMetric(6,12);
             PathMetric(7) + BranchMetric(7,12);
             PathMetric(8) + BranchMetric(8,12);
             PathMetric(9) + BranchMetric(9,12);
             PathMetric(10) + BranchMetric(10,12);
             PathMetric(11) + BranchMetric(11,12);
             PathMetric(12) + BranchMetric(12,12);
             PathMetric(13) + BranchMetric(13,12);
             PathMetric(14) + BranchMetric(14,12);
             PathMetric(15) + BranchMetric(15,12);
             PathMetric(16) + BranchMetric(16,12)];
    [NewPathMetric(12),idx12] = min( cost12 );
    switch idx12
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP12 = [Path,idx12];

    % update P13
    cost13 = [PathMetric(1) + BranchMetric(1,13);
             PathMetric(2) + BranchMetric(2,13);
             PathMetric(3) + BranchMetric(3,13);
             PathMetric(4) + BranchMetric(4,13);
             PathMetric(5) + BranchMetric(5,13);
             PathMetric(6) + BranchMetric(6,13);
             PathMetric(7) + BranchMetric(7,13);
             PathMetric(8) + BranchMetric(8,13);
             PathMetric(9) + BranchMetric(9,13);
             PathMetric(10) + BranchMetric(10,13);
             PathMetric(11) + BranchMetric(11,13);
             PathMetric(12) + BranchMetric(12,13);
             PathMetric(13) + BranchMetric(13,13);
             PathMetric(14) + BranchMetric(14,13);
             PathMetric(15) + BranchMetric(15,13);
             PathMetric(16) + BranchMetric(16,13)];
    [NewPathMetric(13),idx13] = min( cost13 );
    switch idx13
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP13 = [Path,idx13];

    % update P14
    cost14 = [PathMetric(1) + BranchMetric(1,14);
             PathMetric(2) + BranchMetric(2,14);
             PathMetric(3) + BranchMetric(3,14);
             PathMetric(4) + BranchMetric(4,14);
             PathMetric(5) + BranchMetric(5,14);
             PathMetric(6) + BranchMetric(6,14);
             PathMetric(7) + BranchMetric(7,14);
             PathMetric(8) + BranchMetric(8,14);
             PathMetric(9) + BranchMetric(9,14);
             PathMetric(10) + BranchMetric(10,14);
             PathMetric(11) + BranchMetric(11,14);
             PathMetric(12) + BranchMetric(12,14);
             PathMetric(13) + BranchMetric(13,14);
             PathMetric(14) + BranchMetric(14,14);
             PathMetric(15) + BranchMetric(15,14);
             PathMetric(16) + BranchMetric(16,14)];
    [NewPathMetric(14),idx14] = min( cost14 );
    switch idx14
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP14 = [Path,idx14];

    % update P15
    cost15 = [PathMetric(1) + BranchMetric(1,15);
             PathMetric(2) + BranchMetric(2,15);
             PathMetric(3) + BranchMetric(3,15);
             PathMetric(4) + BranchMetric(4,15);
             PathMetric(5) + BranchMetric(5,15);
             PathMetric(6) + BranchMetric(6,15);
             PathMetric(7) + BranchMetric(7,15);
             PathMetric(8) + BranchMetric(8,15);
             PathMetric(9) + BranchMetric(9,15);
             PathMetric(10) + BranchMetric(10,15);
             PathMetric(11) + BranchMetric(11,15);
             PathMetric(12) + BranchMetric(12,15);
             PathMetric(13) + BranchMetric(13,15);
             PathMetric(14) + BranchMetric(14,15);
             PathMetric(15) + BranchMetric(15,15);
             PathMetric(16) + BranchMetric(16,15)];
    [NewPathMetric(15),idx15] = min( cost15 );
    switch idx15
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP15 = [Path,idx15];

    % update P16
    cost16 = [PathMetric(1) + BranchMetric(1,16);
             PathMetric(2) + BranchMetric(2,16);
             PathMetric(3) + BranchMetric(3,16);
             PathMetric(4) + BranchMetric(4,16);
             PathMetric(5) + BranchMetric(5,16);
             PathMetric(6) + BranchMetric(6,16);
             PathMetric(7) + BranchMetric(7,16);
             PathMetric(8) + BranchMetric(8,16);
             PathMetric(9) + BranchMetric(9,16);
             PathMetric(10) + BranchMetric(10,16);
             PathMetric(11) + BranchMetric(11,16);
             PathMetric(12) + BranchMetric(12,16);
             PathMetric(13) + BranchMetric(13,16);
             PathMetric(14) + BranchMetric(14,16);
             PathMetric(15) + BranchMetric(15,16);
             PathMetric(16) + BranchMetric(16,16)];
    [NewPathMetric(16),idx16] = min( cost16 );
    switch idx16
        case 1
            Path = Path1;
        case 2
            Path = Path2;
        case 3
            Path = Path3;
        case 4
            Path = Path4;
        case 5
            Path = Path5;
        case 6
            Path = Path6;
        case 7
            Path = Path7;
        case 8
            Path = Path8;
        case 9
            Path = Path9;
        case 10
            Path = Path10;
        case 11
            Path = Path11;
        case 12
            Path = Path12;
        case 13
            Path = Path13;
        case 14
            Path = Path14;
        case 15
            Path = Path15;
        case 16
            Path = Path16;
    end
    NewP16 = [Path,idx16];

    Path1 = NewP1;
    Path2 = NewP2;
    Path3 = NewP3;
    Path4 = NewP4;
    Path5 = NewP5;
    Path6 = NewP6;
    Path7 = NewP7;
    Path8 = NewP8;
    Path9 = NewP9;
    Path10 = NewP10;
    Path11 = NewP11;
    Path12 = NewP12;
    Path13 = NewP13;
    Path14 = NewP14;
    Path15 = NewP15;
    Path16 = NewP16;
    PathMetric = NewPathMetric;
end
P = [Path1;Path2;Path3;Path4;Path5;Path6;Path7;Path8;...
     Path9;Path10;Path11;Path12;Path13;Path14;Path15;Path16];
% P = [Path1;Path2;Path3;Path4];
[~,indx] = min(PathMetric);
FinalPath = [P(indx,:),indx];
% DecodeBit = zeros(1,numBits*4);

%% Decoding depends on state transitions
% for i = 1 : length(FinalPath)-1
    DecodeBit = reshape(de2bi(FinalPath(1:end-1)-1,4,"left-msb").',1,[]).';
    % switch FinalPath(i)
    %     case 1
    %         DecodeBit( (i-1)*2+1 : i*2) = [0,0];
    %         % DecodeBit(i+1) = 0;
    %     case 2
    %         DecodeBit( (i-1)*2+1 : i*2) = [0,1];
    %     case 3
    %         DecodeBit( (i-1)*2+1 : i*2) = [1,0];
    %     case 4
    %         DecodeBit( (i-1)*2+1 : i*2) = [1,1];
    % end
% end
FinalPath = FinalPath-1;


end