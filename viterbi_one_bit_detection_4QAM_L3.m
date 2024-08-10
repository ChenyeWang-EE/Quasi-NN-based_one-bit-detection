function [DecodeBit,DecodeSym] = viterbi_one_bit_detection_4QAM_L3(zk, channelH,noiseVar)
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
numBits = length(zk);
% % test

for i = 1 : numBits
    NewPathMetric = zeros(16,1);
    % BranchMetric  = distance_branch_one_bit_L2( zk(i),channelH,noiseVar); % distance -log(p)
    BranchMetric = distance_branch_4QAM_L3( zk(i),channelH,noiseVar);
    %% viterbi decode
    % update P1
    cost1 = [PathMetric(1) + BranchMetric(1,1);
             PathMetric(2) + BranchMetric(2,1);
             PathMetric(3) + BranchMetric(3,1);
             PathMetric(4) + BranchMetric(4,1)];
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
    end
    NewP1 = [Path,idx1];
    
    % update P2
    cost2 = [PathMetric(5) + BranchMetric(5,1);
             PathMetric(6) + BranchMetric(6,1);
             PathMetric(7) + BranchMetric(7,1);
             PathMetric(8) + BranchMetric(8,1)];
    [NewPathMetric(2),idx2] = min( cost2 );
    switch idx2
        case 1
            Path = Path5;
        case 2
            Path = Path6;
        case 3
            Path = Path7;
        case 4
            Path = Path8;
    end
    NewP2 = [Path,idx2+4];

    % update P3
    cost3 = [PathMetric(9) + BranchMetric(9,1);
             PathMetric(10) + BranchMetric(10,1);
             PathMetric(11) + BranchMetric(11,1);
             PathMetric(12) + BranchMetric(12,1)];
    [NewPathMetric(3),idx3] = min( cost3 );
    switch idx3
        case 1
            Path = Path9;
        case 2
            Path = Path10;
        case 3
            Path = Path11;
        case 4
            Path = Path12;
    end
    NewP3 = [Path,idx3+8];

    % update P4
    cost4 = [PathMetric(13) + BranchMetric(13,1);
             PathMetric(14) + BranchMetric(14,1);
             PathMetric(15) + BranchMetric(15,1);
             PathMetric(16) + BranchMetric(16,1)];
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
    end
    NewP4 = [Path,idx4+12];

    % update P5
    cost5 = [PathMetric(1) + BranchMetric(1,2);
             PathMetric(2) + BranchMetric(2,2);
             PathMetric(3) + BranchMetric(3,2);
             PathMetric(4) + BranchMetric(4,2)];
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
    end
    NewP5 = [Path,idx5];

    % update P6
    cost6 = [PathMetric(5) + BranchMetric(5,2);
             PathMetric(6) + BranchMetric(6,2);
             PathMetric(7) + BranchMetric(7,2);
             PathMetric(8) + BranchMetric(8,2)];
    [NewPathMetric(6),idx6] = min( cost6 );
    switch idx6
        case 1
            Path = Path5;
        case 2
            Path = Path6;
        case 3
            Path = Path7;
        case 4
            Path = Path8;
    end
    NewP6 = [Path,idx6+4];

    % update P7
    cost7 = [PathMetric(9) + BranchMetric(9,2);
             PathMetric(10) + BranchMetric(10,2);
             PathMetric(11) + BranchMetric(11,2);
             PathMetric(12) + BranchMetric(12,2)];
    [NewPathMetric(7),idx7] = min( cost7 );
    switch idx7
        case 1
            Path = Path9;
        case 2
            Path = Path10;
        case 3
            Path = Path11;
        case 4
            Path = Path12;
    end
    NewP7 = [Path,idx7+8];

    % update P8
    cost8 = [PathMetric(13) + BranchMetric(13,2);
             PathMetric(14) + BranchMetric(14,2);
             PathMetric(15) + BranchMetric(15,2);
             PathMetric(16) + BranchMetric(16,2)];
    [NewPathMetric(8),idx8] = min( cost8 );
    switch idx8
        case 1
            Path = Path13;
        case 2
            Path = Path14;
        case 3
            Path = Path15;
        case 4
            Path = Path16;
    end
    NewP8 = [Path,idx8+12];

    % update P9
    cost9 = [PathMetric(1) + BranchMetric(1,3);
             PathMetric(2) + BranchMetric(2,3);
             PathMetric(3) + BranchMetric(3,3);
             PathMetric(4) + BranchMetric(4,3)];
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
    end
    NewP9 = [Path,idx9];

    % update P10
    cost10 = [PathMetric(5) + BranchMetric(5,3);
              PathMetric(6) + BranchMetric(6,3);
              PathMetric(7) + BranchMetric(7,3);
              PathMetric(8) + BranchMetric(8,3)];
    [NewPathMetric(10),idx10] = min( cost10 );
    switch idx10
        case 1
            Path = Path5;
        case 2
            Path = Path6;
        case 3
            Path = Path7;
        case 4
            Path = Path8;
    end
    NewP10 = [Path,idx10+4];

    % update P11
    cost11 = [PathMetric(9) + BranchMetric(9,3);
              PathMetric(10) + BranchMetric(10,3);
              PathMetric(11) + BranchMetric(11,3);
              PathMetric(12) + BranchMetric(12,3)];
    [NewPathMetric(11),idx11] = min( cost11 );
    switch idx11
        case 1
            Path = Path9;
        case 2
            Path = Path10;
        case 3
            Path = Path11;
        case 4
            Path = Path12;
    end
    NewP11 = [Path,idx11+8];

    % update P12
    cost12 = [PathMetric(13) + BranchMetric(13,3);
              PathMetric(14) + BranchMetric(14,3);
              PathMetric(15) + BranchMetric(15,3);
              PathMetric(16) + BranchMetric(16,3)];
    [NewPathMetric(12),idx12] = min( cost12 );
    switch idx12
        case 1
            Path = Path13;
        case 2
            Path = Path14;
        case 3
            Path = Path15;
        case 4
            Path = Path16;
    end
    NewP12 = [Path,idx12+12];

    % update P13
    cost13 = [PathMetric(1) + BranchMetric(1,4);
              PathMetric(2) + BranchMetric(2,4);
              PathMetric(3) + BranchMetric(3,4);
              PathMetric(4) + BranchMetric(4,4)];
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
    end
    NewP13 = [Path,idx13];

    % update P14
    cost14 = [PathMetric(5) + BranchMetric(5,4);
              PathMetric(6) + BranchMetric(6,4);
              PathMetric(7) + BranchMetric(7,4);
              PathMetric(8) + BranchMetric(8,4)];
    [NewPathMetric(14),idx14] = min( cost14 );
    switch idx14
        case 1
            Path = Path5;
        case 2
            Path = Path6;
        case 3
            Path = Path7;
        case 4
            Path = Path8;
    end
    NewP14 = [Path,idx14+4];

    % update P15
    cost15 = [PathMetric(9) + BranchMetric(9,4);
              PathMetric(10) + BranchMetric(10,4);
              PathMetric(11) + BranchMetric(11,4);
              PathMetric(12) + BranchMetric(12,4)];
    [NewPathMetric(15),idx15] = min( cost15 );
    switch idx15
        case 1
            Path = Path9;
        case 2
            Path = Path10;
        case 3
            Path = Path11;
        case 4
            Path = Path12;
    end
    NewP15 = [Path,idx15+8];

    % update P16
    cost16 = [PathMetric(13) + BranchMetric(13,4);
              PathMetric(14) + BranchMetric(14,4);
              PathMetric(15) + BranchMetric(15,4);
              PathMetric(16) + BranchMetric(16,4)];
    [NewPathMetric(16),idx16] = min( cost16 );
    switch idx16
        case 1
            Path = Path13;
        case 2
            Path = Path14;
        case 3
            Path = Path15;
        case 4
            Path = Path16;
    end
    NewP16 = [Path,idx16+12];

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
% P = [Path1;Path2;Path3;Path4;Path5;Path6;Path7;Path8];
P = [Path1;Path2;Path3;Path4;Path5;Path6;Path7;Path8;...
     Path9;Path10;Path11;Path12;Path13;Path14;Path15;Path16];
[~,indx] = min(PathMetric);
FinalPath = [P(indx,:),indx];
% DecodeBit = zeros(1,numBits*2);

%% Decoding depends on state transitions
% for i = 1 : length(FinalPath)-1
%     switch FinalPath(i)
%         case 1
%             DecodeBit( (i-1)*2+1 : i*2) = [0,0];
%             % DecodeBit(i+1) = 0;
%         case 2
%             DecodeBit( (i-1)*2+1 : i*2) = [0,1];
%         case 3
%             DecodeBit( (i-1)*2+1 : i*2) = [1,0];
%         case 4
%             DecodeBit( (i-1)*2+1 : i*2) = [1,1];
%     end
% end
DecodeSym = ceil(FinalPath/4)-1;
DecodeBit = reshape(de2bi(DecodeSym,2,"left-msb").',1,[]);
end
