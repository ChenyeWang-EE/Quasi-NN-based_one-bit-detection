function DecodeBit  = Viterbi_ISI_SIMO_L2( ykVec,channelH_Vec)
%%% Copyright by Department of Communication Science and Engineering ,Fudan

Path1 = [];
Path2 = [];
Path3 = [];
Path4 = [];
Path5 = [];
Path6 = [];
Path7 = [];
Path8 = [];

PathMetric = zeros(1,8);
PathMetric(1) = 0;
PathMetric(2) = 0;
numBits = length(ykVec);
% branch = branch_2H( chanH );
% test
BranchMetricMtr = distance_branch_SIMO( ykVec,channelH_Vec);
%     branch  = prob_1bitQ(sigma,thetaHat);
for ii = 1 : numBits
    NewPathMetric = zeros(8,1);
    BranchMetric  = BranchMetricMtr(ii,:); % distance -log(p)

    %viterbi decode
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

    %     NewPathMetric
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
for ii = 1:length(FinalPath)-1
    switch FinalPath(ii)
        case 1
            if FinalPath(ii+1) == 7
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
        case 2
            if FinalPath(ii+1) == 3
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
        case 3
            if FinalPath(ii+1) == 1
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
        case 4
            if FinalPath(ii+1) == 5
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
        case 5
            if FinalPath(ii+1) == 3
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
        case 6
            if FinalPath(ii+1) == 7
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
        case 7
            if FinalPath(ii+1) == 5
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
        case 8
            if FinalPath(ii+1) == 1
                DecodeBit(ii) = 0;
            else
                DecodeBit(ii) = 1;
            end
    end
end


end

