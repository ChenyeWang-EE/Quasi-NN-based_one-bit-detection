%% ===========================================
% Function: Equalization forward/backward algorithm
% =================================================
% this algorithm is designed for sequence detection
%  Input : x (length 2^N),lecky<=>prior probability
% Output : result（soft decision）
% priority_p    the prior LLR of x
% noiseVar    the noise variance 噪声方差
%%
function LLR = QNN_detection_SIMO_ISI_Channel(LCk_n, testdata, Re_hMtr,Im_hMtr,flag)
%% upsample = 1;
[len,~] = size(testdata);
Lcky = zeros(len, 1);
% trellisOutput = channelH.*[0 0.7071+0.7071i 0 0.7071-0.7071i;
%     0.7071+0.7071i 0 -0.7071+0.7071i 0;
%     0 -0.7071+0.7071i 0 -0.7071-0.7071i;
%     0.7071-0.7071i 0  -0.7071-0.7071i 0]; %t=T、2T处采样
xkMtr = [ exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(-0.25+0.5)), exp(1j*pi*(0.25+0.5)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)),...
         exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25));
         exp(1j*pi*(0-0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(0.5+0.25)),...
         exp(1j*pi*(1-0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5+0.25)), exp(1j*pi*(-0.5+0.25))];
S = 16;
yprob = zeros(len,S);
for ii = 1 : length(testdata)
    % [yprob(ii,:)] = QNN_SIMO_ISI_Channel(testdata(ii,:), Re_hMtr,Im_hMtr,xkMtr);
    [yprob(ii,:)] = QNN_SIMO_ISI_Channel_2(testdata(ii,:), Re_hMtr,Im_hMtr,xkMtr);
end
yprob = yprob + eps;
gammakpost = zeros(8,8,len);
M = [0 0 0 0 0 0 1 2;
     0 0 3 4 0 0 0 0;
     5 6 0 0 0 0 0 0;
     0 0 0 0 7 8 0 0;
     0 0 9 10 0 0 0 0;
     0 0 0 0 0 0 11 12;
     0 0 0 0 13 14 0 0;
     15 16 0 0 0 0 0 0];
R = (1:1:16);
[~,pos]=ismember(R,M);
A = zeros(size(M));
for mm = 1:len
    A(pos)=yprob(mm,:);
    gammakpost(:,:,mm) = A;
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

noZeroPos = find(matrices~=0);
PkPost = zeros(8,8,len);
bk1 =  ones(8,1,len+1);  %bk

gammak = zeros(8,8,len);
for mm =1:len
    %         a = exp(Apos.*LCk(mm))/(1 + exp(LCk(mm))).*matrices;
    b = exp(Apos.*LCk_n(mm))/(1 + exp(LCk_n(mm))).*matrices;
    c = b(noZeroPos)./1;
    c = c/sum(c(1:2));
    d = gammak(:,:,mm);
    d(noZeroPos) = c;
    gammak(:,:,mm) = d;
end



for ii = len:-1:1
    PkPost(:,:,ii) = gammakpost(:,:,ii).*gammak(:,:,ii); % Pk,gamma构成的P矩阵
    bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
    bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));

end
fk = ones(8,1); %f0
%     fk = [1;0;0;0];
%     Pkinit = rand(4,4).*matrices;
%     fk = Pkinit'*fk; %f1
%     Lcky(1) = log(fk'*(Apos.*Pkinit)*bk1(:,:,2)/(fk'*(Aneg.*Pkinit)*bk1(:,:,2)));
Pk1 = gammakpost(:,:,1).*gammak(:,:,1);
Lcky(1) = log(fk'*(Apos.*Pk1)*bk1(:,:,2)/(fk'*(Aneg.*Pk1)*bk1(:,:,2)));

for kk = 2 : len

    PkPrior = gammakpost(:,:,(kk-1)).*gammak(:,:,kk-1);  % P,gamma构成的P矩阵
    fk = PkPrior'*fk;   % 前向更新 fk
    fk = fk / sum(fk);
    Pk=gammakpost(:,:,kk).*gammak(:,:,kk);
    Lcky(kk) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));

end
if numel(Lcky(isnan(Lcky)))>0
    disp('NaN')
    keyboard
end
if flag == 0
    Lcky = Lcky(2:end-1);
    LLR = Lcky;
    LLR(LLR>0)=1;
    LLR(LLR<0)=0;
else
    LLR = Lcky;
end
end