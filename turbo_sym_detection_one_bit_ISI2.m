%% ===========================================
% Function: Equalization forward/backward algorithm
% =================================================
% this algorithm is designed for sequence detection
%  Input : x (length 2^N),lecky<=>prior probability
% Output : result（soft decision）
% priority_p    the prior LLR of x
% noiseVar    the noise variance 噪声方差
%%
function LLR = turbo_sym_detection_one_bit_ISI2(zk, Leckp, noiseVar, channelH,flag)
%% upsample = 1;
len = length(zk);
Lcky = zeros(len, 1);
yprob = zeros(8,8,len);
    trellisOutput = statechangeoutput(channelH);

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
% noiseVar2 = noiseVar;
%     sqrtnoiseVar = pi*noiseVar2;
gammakpost = zeros(8,8,len);
PkPost = zeros(8,8,len);
bk1 =  ones(8,1,len+1);  %bk
sigma = sqrt(noiseVar);
%yprob
for ii = 1 : len
    switch zk(ii)
        case (1+1i) ./ sqrt(2)
            yprob(:,:,ii) = qfunc(-real(trellisOutput)./sigma).*qfunc(-imag(trellisOutput)./sigma);
        case (1-1i) ./ sqrt(2) 
            yprob(:,:,ii) = qfunc(-real(trellisOutput)./sigma).*( 1-qfunc(-imag(trellisOutput)./sigma) );
        case (-1-1i) ./ sqrt(2)
            yprob(:,:,ii) = ( 1-qfunc(-real(trellisOutput)./sigma) ).*( 1-qfunc(-imag(trellisOutput)./sigma) );
        case (-1+1i) ./ sqrt(2)
            yprob(:,:,ii) = ( 1-qfunc(-real(trellisOutput)./sigma) ).*qfunc(-imag(trellisOutput)./sigma);
    end
end
yprob = yprob + eps;
for ii = len:-1:1
    gammakpost(:,:,ii) = exp(Apos.*Leckp(ii))/(1 + exp(Leckp(ii))).*matrices; % P(xk=xi,j)
%     switch zk(ii)
%         case 1+1i
%             yprob = qfunc(-real(trellisOutput)./sigma).*qfunc(-imag(trellisOutput)./sigma);
%         case 1-1i
%             yprob = qfunc(-real(trellisOutput)./sigma).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1-1i
%             yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1+1i
%             yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*qfunc(-imag(trellisOutput)./sigma);
%     end
    PkPost(:,:,ii) = gammakpost(:,:,ii).* yprob(:,:,ii); % Pk,gamma构成的P矩阵

    %     PkPost(:,:,ii) = gammakpost(:,:,ii).* (exp(-(abs(zk(ii)-trellisOutput)).^2/noiseVar2)); % Pk,gamma构成的P矩阵

    bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
    bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));

end
    fk = ones(8,1); %f0
% fk = [1;0;0;0];
%     Pkinit = rand(4,4).*matrices;
%     fk = Pkinit'*fk; %f1
%     Lcky(1) = log(fk'*(Apos.*Pkinit)*bk1(:,:,2)/(fk'*(Aneg.*Pkinit)*bk1(:,:,2)));
gammak = exp(Apos.*Leckp(1))/(1 + exp(Leckp(1))).*matrices;
% switch zk(1)
%     case 1+1i
%         yprob = qfunc(-real(trellisOutput)./sigma).*qfunc(-imag(trellisOutput)./sigma);
%     case 1-1i
%         yprob = qfunc(-real(trellisOutput)./sigma).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%     case -1-1i
%         yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%     case -1+1i
%         yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*qfunc(-imag(trellisOutput)./sigma);
% end
% Pk1=gammak.*exp(-(abs(zk(1)-trellisOutput)).^2/noiseVar2);
Pk1=gammak.* yprob(:,:,1);
tmp = log(fk'*(Apos.*Pk1)*bk1(:,:,2)/(fk'*(Aneg.*Pk1)*bk1(:,:,2)));
Lcky(1)=tmp;
for kk = 2 : len

    gammakprior = exp(Apos.*Leckp(kk-1))/(1 + exp(Leckp(kk-1))).*matrices; % P(xk=xi,j)
%     switch zk(kk-1)
%         case 1+1i
%             yprob = qfunc(-real(trellisOutput)./sigma).*qfunc(-imag(trellisOutput)./sigma);
%         case 1-1i
%             yprob = qfunc(-real(trellisOutput)./sigma).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1-1i
%             yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1+1i
%             yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*qfunc(-imag(trellisOutput)./sigma);
%     end
        PkPrior = gammakprior.*yprob(:,:,kk-1);  % P,gamma构成的P矩阵

%     PkPrior = gammakprior.*(exp(-(abs(zk(kk-1)-trellisOutput)).^2/noiseVar2));  % P,gamma构成的P矩阵
    fk = PkPrior'*fk;   % 前向更新 fk
    fk = fk / sum(fk);
    gammak = exp(Apos.*Leckp(kk))/(1 + exp(Leckp(kk))).*matrices;
%     switch zk(kk)
%         case 1+1i
%             yprob = qfunc(-real(trellisOutput)./sigma).*qfunc(-imag(trellisOutput)./sigma);
%         case 1-1i
%             yprob = qfunc(-real(trellisOutput)./sigma).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1-1i
%             yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1+1i
%             yprob = ( 1-qfunc(-real(trellisOutput)./sigma) ).*qfunc(-imag(trellisOutput)./sigma);
%     end
    %     Pk=gammak.*exp(-(abs(zk(kk)-trellisOutput)).^2/noiseVar2);
    Pk=gammak.*yprob(:,:,kk);
    Lcky(kk) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));

end

%     Lcky = Lcky(3:end-2);
if numel(Lcky(isnan(Lcky)))>0
    disp('NaN')
    keyboard
end

if flag == 1
    LLR = Lcky;
else
    Lcky = Lcky(2:end-1);
    LLR = Lcky;
    LLR(LLR>0)=1;
    LLR(LLR<0)=0;
end
end

