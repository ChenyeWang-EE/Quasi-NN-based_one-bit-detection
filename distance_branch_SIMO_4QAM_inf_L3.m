function BranchMetricMtr = distance_branch_SIMO_4QAM_inf_L3( ykVec, channelH_Vec)

% len = length(branch);
[antennaNum,len] = size(ykVec);
branchTmp = zeros(16,4,len.*antennaNum);
branch = zeros(16,4,len);

for ii = 1 : antennaNum
    % trellisOutput(8*ii-7 : 8*ii,:) = statechangeoutput(channelH_Vec(ii,:)); %t=T、2T处采样
    trellisOutput(16*ii-15 : 16*ii,:) = state_tran_output_4QAM_L3(channelH_Vec(ii,:)); %t=T、2T处采样
end
% distance = zeros(1,len);
% for kk = 1:len
%     % distance(kk) = abs(branch(kk)).^2-2*real(sig'*branch(kk));
%     distance(kk) = norm(branch(kk)-sig,2).^2;
%     % distance(kk) = ( abs(branch(kk)-sig) ).^2/noiseVar;
%     % distance3(kk) = ( branch(kk) - sig )'*( branch(kk) - sig );
% end

for ii = 1 : len
    for jj = 1 : antennaNum
        getnormMtr = ykVec(jj,ii) - trellisOutput(16*jj-15 : 16*jj,:);
        branchTmp(:,:,(ii-1).*antennaNum+jj) = getnormMtr.*conj(getnormMtr) + eps; % Pk,gamma构成的P矩阵
        branch(:,:,ii) =  branchTmp(:,:,(ii-1).*antennaNum+jj) + branch(:,:,ii);
    end
    % tmp = branch(:,:,ii);
    % BranchMetricMtr(ii,:) = tmp(pos);
end
BranchMetricMtr = branch;
end
