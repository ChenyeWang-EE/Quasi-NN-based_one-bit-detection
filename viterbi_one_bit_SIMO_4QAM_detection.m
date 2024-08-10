function DecodeBit = viterbi_one_bit_SIMO_4QAM_detection(zkVec, channelH_Vec,noiseVar)
%% upsample = 1;
% len = length(zk);
% yprob = zeros(8,8,len);

Path1 = [];
Path2 = [];
Path3 = [];
Path4 = [];

PathMetric = zeros(1,4);
% PathMetric(1) = 0;
% PathMetric(2) = 0;
% numBits = length(zkVec);
numBits = size(zkVec,2);
% % test
BranchMetricMtr = distance_branch_SIMO_4QAM( zkVec,channelH_Vec,noiseVar);
for i = 1 : numBits
    NewPathMetric = zeros(4,1);
    % BranchMetric  = distance_branch_one_bit_L2( zk(i),channelH,noiseVar); % distance -log(p)
    % BranchMetric = distance_branch_4QAM( zk,channelH,noiseVar);
    BranchMetric = BranchMetricMtr(:,:,i); 
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
    cost2 = [PathMetric(1) + BranchMetric(1,2);
             PathMetric(2) + BranchMetric(2,2);
             PathMetric(3) + BranchMetric(3,2);
             PathMetric(4) + BranchMetric(4,2)];
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
    end
    NewP2 = [Path,idx2];

    % update P3
    cost3 = [PathMetric(1) + BranchMetric(1,3);
             PathMetric(2) + BranchMetric(2,3);
             PathMetric(3) + BranchMetric(3,3);
             PathMetric(4) + BranchMetric(4,3)];
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
    end
    NewP3 = [Path,idx3];

    % update P4
    cost4 = [PathMetric(1) + BranchMetric(1,4);
             PathMetric(2) + BranchMetric(2,4);
             PathMetric(3) + BranchMetric(3,4);
             PathMetric(4) + BranchMetric(4,4)];
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
    NewP4 = [Path,idx4];

    Path1 = NewP1;
    Path2 = NewP2;
    Path3 = NewP3;
    Path4 = NewP4;
    % Path5 = NewP5;
    % Path6 = NewP6;
    % Path7 = NewP7;
    % Path8 = NewP8;
    PathMetric = NewPathMetric;
end
% P = [Path1;Path2;Path3;Path4;Path5;Path6;Path7;Path8];
P = [Path1;Path2;Path3;Path4];
[~,indx] = min(PathMetric);
FinalPath = [P(indx,:),indx];
DecodeBit = zeros(1,numBits*2);

%% Decoding depends on state transitions
for i = 1 : length(FinalPath)-1
    switch FinalPath(i)
        case 1
            DecodeBit( (i-1)*2+1 : i*2) = [0,0];
            % DecodeBit(i+1) = 0;
        case 2
            DecodeBit( (i-1)*2+1 : i*2) = [0,1];
        case 3
            DecodeBit( (i-1)*2+1 : i*2) = [1,0];
        case 4
            DecodeBit( (i-1)*2+1 : i*2) = [1,1];
        % % case 5
        % %     if FinalPath(i+1) == 3
        % %         DecodeBit(i) = 0;
        % %     else
        % %         DecodeBit(i) = 1;
        % %     end
        % % case 6
        % %     if FinalPath(i+1) == 7
        % %         DecodeBit(i) = 0;
        % %     else
        % %         DecodeBit(i) = 1;
        % %     end
        % % case 7
        % %     if FinalPath(i+1) == 5
        % %         DecodeBit(i) = 0;
        % %     else
        % %         DecodeBit(i) = 1;
        % %     end
        % % case 8
        % %     if FinalPath(i+1) == 1
        % %         DecodeBit(i) = 0;
        % %     else
        % %         DecodeBit(i) = 1;
        % %     end
    end
end



end
