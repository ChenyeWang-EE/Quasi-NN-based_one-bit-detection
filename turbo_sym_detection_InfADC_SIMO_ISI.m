%% ===========================================
% Function: Equalization forward/backward algorithm
% =================================================
% this algorithm is designed for sequence detection
%  Input : x (length 2^N),lecky<=>prior probability
% Output : result（soft decision）
% priority_p    the prior LLR of x
% noiseVar    the noise variance 噪声方差
%%
function LLR = turbo_sym_detection_InfADC_SIMO_ISI(ykVec, Leckp, noiseVar, channelH_Vec,flag)
%% upsample = 1;
% len = length(zk);
[antennaNum,len] = size(ykVec);
Lcky = zeros(len, 1);
yprobTmp = zeros(8,8,len.*antennaNum);
% yprob = ones(8,8,len);
yprob = zeros(8,8,len);
for ii = 1 : antennaNum
    trellisOutput(8*ii-7 : 8*ii,:) = statechangeoutput(channelH_Vec(ii,:)); %t=T、2T处采样
end

     Apos = [0 0 0 0 0 0 0 1;
             0 0 0 1 0 0 0 0;
             0 1 0 0 0 0 0 0;
             0 0 0 0 0 1 0 0;
             0 0 0 1 0 0 0 0;
             0 0 0 0 0 0 0 1;
             0 0 0 0 0 1 0 0;
             0 1 0 0 0 0 0 0]; %  A(+1)
     Aneg = [0 0 0 0 0 0 1 0;
             0 0 1 0 0 0 0 0;
             1 0 0 0 0 0 0 0;
             0 0 0 0 1 0 0 0;
             0 0 1 0 0 0 0 0;
             0 0 0 0 0 0 1 0;
             0 0 0 0 1 0 0 0;
             1 0 0 0 0 0 0 0];  %  A(-1)
 matrices = [0 0 0 0 0 0 1 1;
             0 0 1 1 0 0 0 0;
             1 1 0 0 0 0 0 0;
             0 0 0 0 1 1 0 0;
             0 0 1 1 0 0 0 0;
             0 0 0 0 0 0 1 1;
             0 0 0 0 1 1 0 0;
             1 1 0 0 0 0 0 0];     % for the Pk


%     noiseVar2 = noiseVar+db2pow(-1); %干扰功率加噪声功率
noiseVar2 = noiseVar;
%     sqrtnoiseVar = pi*noiseVar2;
gammakpost = zeros(8,8,len);
PkPost = zeros(8,8,len);
bk1 =  ones(8,1,len+1);  %bk
% sigma = sqrt(noiseVar);
%yprob
for ii = 1 : len
    for jj = 1 : antennaNum
        % switch zkVec(jj,ii)
        %     case (1+1i) ./ sqrt(2)
        %         yprobTmp(:,:,(ii-1).*antennaNum+jj) = qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma).*qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma);
        %     case (1-1i) ./ sqrt(2)
        %         yprobTmp(:,:,(ii-1).*antennaNum+jj) = qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma).*( 1-qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma) );
        %     case (-1-1i) ./ sqrt(2)
        %         yprobTmp(:,:,(ii-1).*antennaNum+jj) = ( 1-qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma) ).*( 1-qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma) );
        %     case (-1+1i) ./ sqrt(2)
        %         yprobTmp(:,:,(ii-1).*antennaNum+jj) = ( 1-qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma) ).*qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma);
        % end
        yprobTmp(:,:,(ii-1).*antennaNum+jj) = exp(-(abs(ykVec(jj,ii)-trellisOutput(8*jj-7 : 8*jj,:))).^2/noiseVar2) + eps; % Pk,gamma构成的P矩阵
        % yprob(:,:,ii) = yprobTmp(:,:,(ii-1).*antennaNum+jj).*yprob(:,:,ii);
        yprob(:,:,ii) = log( yprobTmp(:,:,(ii-1).*antennaNum+jj) ) + yprob(:,:,ii);
    end
end
yprob = exp(yprob);
for ii = len:-1:1
    gammakpost(:,:,ii) = exp(Apos.*Leckp(ii))/(1 + exp(Leckp(ii))).*matrices; % P(xk=xi,j)

    PkPost(:,:,ii) = gammakpost(:,:,ii).* yprob(:,:,ii); % Pk,gamma构成的P矩阵

    % PkPost(:,:,ii) = gammakpost(:,:,ii).* (exp(-(abs(ykVec(ii)-trellisOutput)).^2/noiseVar2)); % Pk,gamma构成的P矩阵

    bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
    bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));

end
    fk = ones(8,1); %f0
% fk = [1;0;0;0];
%     Pkinit = rand(4,4).*matrices;
%     fk = Pkinit'*fk; %f1
%     Lcky(1) = log(fk'*(Apos.*Pkinit)*bk1(:,:,2)/(fk'*(Aneg.*Pkinit)*bk1(:,:,2)));
gammak = exp(Apos.*Leckp(1))/(1 + exp(Leckp(1))).*matrices;

% Pk1=gammak.*exp(-(abs(zk(1)-trellisOutput)).^2/noiseVar2);
Pk1=gammak.* yprob(:,:,1);
tmp = log(fk'*(Apos.*Pk1)*bk1(:,:,2)/(fk'*(Aneg.*Pk1)*bk1(:,:,2)));
Lcky(1)=tmp;
for kk = 2 : len

    gammakprior = exp(Apos.*Leckp(kk-1))/(1 + exp(Leckp(kk-1))).*matrices; % P(xk=xi,j)

        PkPrior = gammakprior.*yprob(:,:,kk-1);  % P,gamma构成的P矩阵

%     PkPrior = gammakprior.*(exp(-(abs(zk(kk-1)-trellisOutput)).^2/noiseVar2));  % P,gamma构成的P矩阵
    fk = PkPrior'*fk;   % 前向更新 fk
    fk = fk / sum(fk);
    gammak = exp(Apos.*Leckp(kk))/(1 + exp(Leckp(kk))).*matrices;

    %     Pk=gammak.*exp(-(abs(zk(kk)-trellisOutput)).^2/noiseVar2);
    Pk=gammak.*yprob(:,:,kk);
    Lcky(kk) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));

end

if flag == 1
    % turbo
    LLR = Lcky;
else
    % % GFSK+ISI
    Lcky = Lcky(2:end-1);
    LLR = Lcky;
    LLR(LLR>0)=1;
    LLR(LLR<0)=0;
end
end

