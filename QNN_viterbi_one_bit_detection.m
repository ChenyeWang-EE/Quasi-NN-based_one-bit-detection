function DecodeBit = QNN_viterbi_one_bit_detection(zk,Re_h,Im_h)
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

PathMetric = zeros(1,8);
% PathMetric(1) = 0;
% PathMetric(2) = 0;
numBits = length(zk);
% compute branchMetric using QNN
distanceMtr = QNN_distance_branch_one_bit( zk,Re_h,Im_h); % numBits × 16 Matrix

for i = 1:numBits
    NewPathMetric = zeros(8,1);
    % compute branchMetric using QNN
    BranchMetric  = distanceMtr(i,:); % distance -log(p)

    %% viterbi decoding
    % update P1
    if (PathMetric(3) + BranchMetric(5)) <= (PathMetric(8) + BranchMetric(15))
        NewPathMetric(1) = PathMetric(3) + BranchMetric(5);
        NewP1 = [Path3,3];
    else
        NewPathMetric(1) = PathMetric(8) + BranchMetric(15);
        NewP1 = [Path8,8];
    end

    % update P2
    if (PathMetric(3) + BranchMetric(6)) <= (PathMetric(8) + BranchMetric(16))
        NewPathMetric(2) = PathMetric(3) + BranchMetric(6);
        NewP2 = [Path3,3];
    else
        NewPathMetric(2) = PathMetric(8) + BranchMetric(16);
        NewP2 = [Path8,8];
    end

    % update P3
    if (PathMetric(2) + BranchMetric(3)) <= (PathMetric(5) + BranchMetric(9))
        NewPathMetric(3) = PathMetric(2) + BranchMetric(3);
        NewP3 = [Path2,2];
    else
        NewPathMetric(3) = PathMetric(5) + BranchMetric(9);
        NewP3 = [Path5,5];
    end

    % update P4
    if (PathMetric(2) + BranchMetric(4)) <= (PathMetric(5) + BranchMetric(10))
        NewPathMetric(4) = PathMetric(2) + BranchMetric(4);
        NewP4 = [Path2,2];
    else
        NewPathMetric(4) = PathMetric(5) + BranchMetric(10);
        NewP4 = [Path5,5];
    end

    % update P5
    if (PathMetric(4) + BranchMetric(7)) <= (PathMetric(7) + BranchMetric(13))
        NewPathMetric(5) = PathMetric(4) + BranchMetric(7);
        NewP5 = [Path4,4];
    else
        NewPathMetric(5) = PathMetric(7) + BranchMetric(13);
        NewP5 = [Path7,7];
    end

    % update P6
    if (PathMetric(4) + BranchMetric(8)) <= (PathMetric(7) + BranchMetric(14))
        NewPathMetric(6) = PathMetric(4) + BranchMetric(8);
        NewP6 = [Path4,4];
    else
        NewPathMetric(6) = PathMetric(7) + BranchMetric(14);
        NewP6 = [Path7,7];
    end

    % update P7
    if (PathMetric(1) + BranchMetric(1)) <= (PathMetric(6) + BranchMetric(11))
        NewPathMetric(7) = PathMetric(1) + BranchMetric(1);
        NewP7 = [Path1,1];
    else
        NewPathMetric(7) = PathMetric(6) + BranchMetric(11);
        NewP7 = [Path6,6];
    end

    % update P8
    if (PathMetric(1) + BranchMetric(2)) <= (PathMetric(6) + BranchMetric(12))
        NewPathMetric(8) = PathMetric(1) + BranchMetric(2);
        NewP8 = [Path1,1];
    else
        NewPathMetric(8) = PathMetric(6) + BranchMetric(12);
        NewP8 = [Path6,6];
    end

    Path1 = NewP1;
    Path2 = NewP2;
    Path3 = NewP3;
    Path4 = NewP4;
    Path5 = NewP5;
    Path6 = NewP6;
    Path7 = NewP7;
    Path8 = NewP8;
    PathMetric = NewPathMetric;
end
P = [Path1;Path2;Path3;Path4;Path5;Path6;Path7;Path8];
[~,indx] = min(PathMetric);
FinalPath = [P(indx,:),indx];
DecodeBit = zeros(1,numBits);

%% Decoding depends on state transitions
for i = 1:length(FinalPath)-1
    switch FinalPath(i)
        case 1
            if FinalPath(i+1) == 7
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
        case 2
            if FinalPath(i+1) == 3
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
        case 3
            if FinalPath(i+1) == 1
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
        case 4
            if FinalPath(i+1) == 5
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
        case 5
            if FinalPath(i+1) == 3
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
        case 6
            if FinalPath(i+1) == 7
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
        case 7
            if FinalPath(i+1) == 5
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
        case 8
            if FinalPath(i+1) == 1
                DecodeBit(i) = 0;
            else
                DecodeBit(i) = 1;
            end
    end
end


end
