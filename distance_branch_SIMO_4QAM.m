function [ BranchMetricMtr ] = distance_branch_SIMO_4QAM( zkVec,channelH_Vec,noiseVar)


[antennaNum,len] = size(zkVec);
yprobTmp = zeros(4,4,len.*antennaNum);
yprob = zeros(4,4,len);
BranchMetricMtr = zeros(4,4,len);
for ii = 1 : antennaNum
    trellisOutput(4*ii-3 : 4*ii,:) = state_tran_output_4QAM(channelH_Vec(ii,:)); %t=T、2T处采样
end
% trellisOutput = statechangeoutput(channelH);
% trellisOutput = state_tran_output_4QAM(channelH_Vec);
sigma = sqrt(noiseVar);
%yprob
for ii = 1 : len
    for jj = 1 : antennaNum
        switch zkVec(jj,ii)
            case (1+1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = qfunc(-real(trellisOutput(4*jj-3 : 4*jj,:))./sigma).*qfunc(-imag(trellisOutput(4*jj-3 : 4*jj,:))./sigma);
            case (1-1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = qfunc(-real(trellisOutput(4*jj-3 : 4*jj,:))./sigma).*( 1-qfunc(-imag(trellisOutput(4*jj-3 : 4*jj,:))./sigma) );
            case (-1-1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = ( 1-qfunc(-real(trellisOutput(4*jj-3 : 4*jj,:))./sigma) ).*( 1-qfunc(-imag(trellisOutput(4*jj-3 : 4*jj,:))./sigma) );
            case (-1+1i) ./ sqrt(2)
                yprobTmp(:,:,(ii-1).*antennaNum+jj) = ( 1-qfunc(-real(trellisOutput(4*jj-3 : 4*jj,:))./sigma) ).*qfunc(-imag(trellisOutput(4*jj-3 : 4*jj,:))./sigma);
        end
        yprobTmp(:,:,(ii-1).*antennaNum+jj) = yprobTmp(:,:,(ii-1).*antennaNum+jj) + eps;
        % yprob(:,:,ii) = yprobTmp(:,:,(ii-1).*antennaNum+jj).*yprob(:,:,ii);
        yprob(:,:,ii) = log( yprobTmp(:,:,(ii-1).*antennaNum+jj) ) + yprob(:,:,ii) ;
    end
    % tmp = yprob(:,:,ii);
    BranchMetricMtr(:,:,ii) = - yprob(:,:,ii) ;
end
% yprob = yprob + eps;
% % S = [0 0 0 0 0 0 1 2;
% %     0 0 3 4 0 0 0 0;
% %     5 6 0 0 0 0 0 0;
% %     0 0 0 0 7 8 0 0;
% %     0 0 9 10 0 0 0 0;
% %     0 0 0 0 0 0 11 12;
% %     0 0 0 0 13 14 0 0;
% %     15 16 0 0 0 0 0 0];
% % R = (1:1:16);
% % [~,pos]=ismember(R,S);
% % branch = yprob(pos);
% distance = zeros(1,len);
% for kk = 1:len
%     % distance(kk) = abs(branch(kk)).^2-2*real(sig'*branch(kk));
%     distance(kk) = norm(branch(kk)-sig,2).^2;
%     % distance(kk) = ( abs(branch(kk)-sig) ).^2/noiseVar;
%     % distance3(kk) = ( branch(kk) - sig )'*( branch(kk) - sig );
%     distance(kk)
% end
% distance = -log(yprob);

end
