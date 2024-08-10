function [Re_hMtr,Im_hMtr, lossEpoch] = forward_update_SIMO_Channel(ResigVec, ImsigVec, labelVec, beta, lr, batch_size,testBits)

transition = [0 0.7071+0.7071i 0 0.7071-0.7071i;
              0.7071+0.7071i 0 -0.7071+0.7071i 0;
              0 -0.7071+0.7071i 0 -0.7071-0.7071i;
              0.7071-0.7071i 0  -0.7071-0.7071i 0];
S = [0 2 0 8;
     5 0 3 0;
     0 6 0 4;
     1 0 7 0];
R = (1:1:8);
[~,pos]=ismember(R,S);
xkVec = transition(pos);
% M = 8;
[~,N] = size(ResigVec);
loss = zeros(1,batch_size);
dLdReh = zeros(N,batch_size);
dLdImh = zeros(N,batch_size);
Re_hMtr = zeros(N,1);
Im_hMtr = zeros(N,1);
% Re_hMtr = randn(N,1)+randn(N,1).*1i;
% Im_hMtr = randn(N,1)+randn(N,1).*1i;
% Re_v = -[0.7;0.7;-0.7;-0.7;0.7;-0.7;-0.7;0.7]/0.3162;
% Im_v = -[-0.7;0.7;0.7;-0.7;0.7;0.7;-0.7;-0.7]/0.3162;
iternum = floor(testBits / batch_size);
lossAvg = zeros(1,iternum);
epoch = 50;
vtRe = zeros(N,1); %initial gradient
vtIm = zeros(N,1);
lossEpoch = zeros(1,epoch);
for idx = 1 : epoch
    for jj = 1 : iternum
        % sumdLdReh = zeros(N,1);
        % sumdLdImh = zeros(N,1);
        for ii = 1 : batch_size
            [dLdReh(:,ii),dLdImh(:,ii),loss(ii)] = GenGradient_SIMO_Channel(ResigVec( ii + (jj-1)*batch_size , : ),ImsigVec( ii + (jj-1)*batch_size , : ),Re_hMtr,Im_hMtr,labelVec( ii + (jj-1)*batch_size),xkVec.');
            % sumdLdReh = sumdLdReh + dLdReh(:,ii);
            % sumdLdImh = sumdLdImh + dLdImh(:,ii);
        end
        % gradRe = sumdLdReh ./ batch_size;
        % gradIm = sumdLdImh ./ batch_size;
        gradRe = mean(dLdReh,2);
        gradIm = mean(dLdImh,2);
        lossAvg(jj) = mean(loss);
        % Gradient Descent with Momentum
        vtRe = beta.* vtRe + (1 - beta).*(gradRe);
        vtIm = beta.* vtIm + (1 - beta).*(gradIm);

%         if jj>10 && ( abs(lossAvg(jj)-lossAvg(jj-10)) / abs(lossAvg(jj)) < 1e-5 )
%             break;
%         end
%         Re_vMtr = Re_vMtr - lr.*gradRe;
%         Im_vMtr = Im_vMtr - lr.*gradIm;
        Re_hMtr = Re_hMtr - lr.*vtRe;
        Im_hMtr = Im_hMtr - lr.*vtIm;
    end
    lossEpoch(idx) = mean(lossAvg);
end

end