function BranchMetricMtr = distance_branch_one_bit_SIMO_L2( zkVec,channelH_Vec,noiseVar)

[antennaNum,len] = size(zkVec);
yprobTmp = zeros(8,8,len.*antennaNum);
yprob = zeros(8,8,len);
for ii = 1 : antennaNum
    trellisOutput(8*ii-7 : 8*ii,:) = statechangeoutput(channelH_Vec(ii,:)); %t=T、2T处采样
end
sigma = sqrt(noiseVar);
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
%yprob
for ii = 1 : len
    for jj = 1 : antennaNum
        switch zkVec(jj,ii)
            case (1+1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma).*qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma);
            case (1-1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma).*( 1-qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma) );
            case (-1-1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = ( 1-qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma) ).*( 1-qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma) );
            case (-1+1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = ( 1-qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma) ).*qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma);
        end
        yprobTmp(:,:,(ii-1).*antennaNum+jj) = yprobTmp(:,:,(ii-1).*antennaNum+jj) + eps;
        % yprob(:,:,ii) = yprobTmp(:,:,(ii-1).*antennaNum+jj).*yprob(:,:,ii);
        yprob(:,:,ii) = log( yprobTmp(:,:,(ii-1).*antennaNum+jj) ) + yprob(:,:,ii) ;
    end
    tmp = yprob(:,:,ii);
    BranchMetricMtr(ii,:) = -tmp(pos);
end
% yprob = exp(yprob);


% distance = zeros(1,len);
% for kk = 1:len
%     % distance(kk) = abs(branch(kk)).^2-2*real(sig'*branch(kk));
%     distance(kk) = norm(branch(kk)-sig,2).^2;
%     % distance(kk) = ( abs(branch(kk)-sig) ).^2/noiseVar;
%     % distance3(kk) = ( branch(kk) - sig )'*( branch(kk) - sig );
%     distance(kk)
% end
% distance = branch;

end
