function BranchMetricMtr = distance_branch_SIMO( ykVec,channelH_Vec)

% len = length(branch);
[antennaNum,len] = size(ykVec);
branchTmp = zeros(8,8,len.*antennaNum);
branch = zeros(8,8,len);
S = [0 0 0 0 0 0 1 2;
    0 0 3 4 0 0 0 0;
    5 6 0 0 0 0 0 0;
    0 0 0 0 7 8 0 0;
    0 0 9 10 0 0 0 0;
    0 0 0 0 0 0 11 12;
    0 0 0 0 13 14 0 0;
    15 16 0 0 0 0 0 0];
R = (1:1:16);
[~,pos]=ismember(R,S);
BranchMetricMtr = zeros(len,16);
for ii = 1 : antennaNum
    trellisOutput(8*ii-7 : 8*ii,:) = statechangeoutput(channelH_Vec(ii,:)); %t=T、2T处采样
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
        getnormMtr = ykVec(jj,ii) - trellisOutput(8*jj-7 : 8*jj,:);
        branchTmp(:,:,(ii-1).*antennaNum+jj) = getnormMtr.*conj(getnormMtr) + eps; % Pk,gamma构成的P矩阵
        branch(:,:,ii) =  branchTmp(:,:,(ii-1).*antennaNum+jj) + branch(:,:,ii);
    end
    tmp = branch(:,:,ii);
    BranchMetricMtr(ii,:) = tmp(pos);
end

end
