function [Re_vMtr,Im_vMtr, loss] = forward_update_SIMO_ISI(ResigVec, ImsigVec, labelVec, beta, lr, batch_size,testBits)

M = 16;
[~,N] = size(ResigVec);
loss = zeros(1,batch_size);
% loss3 = zeros(1,batch_size);
dLdRev = zeros(M,N,batch_size);
dLdImv = zeros(M,N,batch_size);
Re_vMtr = zeros(M,N);
Im_vMtr = zeros(M,N);
% Re_v = -[0.7;0.7;-0.7;-0.7;0.7;-0.7;-0.7;0.7]/0.3162;
% Im_v = -[-0.7;0.7;0.7;-0.7;0.7;0.7;-0.7;-0.7]/0.3162;
iternum = floor((testBits-1) / batch_size);
lossAvg = zeros(1,iternum);
% lossAvg2 = zeros(1,iternum);
epoch = 50;
loss2 = zeros(1,epoch);
% loss4 = zeros(1,epoch);
vtRe = zeros(M,N); %initial gradient
vtIm = zeros(M,N);
% %%---------------------%%
% labelyj = zeros(testBits-1,16);
% yprobTmp = zeros(8,8,(testBits-1).*N);
% yprob = ones(8,8,testBits-1);
% for ii = 1 : N
%     trellisOutput(8*ii-7 : 8*ii,:) = statechangeoutput(channelH_Vec(ii,:)); %t=T、2T处采样
% end
% %yprob
% zkVec = (ResigVec + 1i*ImsigVec).';
% sigma = sqrt(noiseVar);
% 
% S = [0 0 0 0 0 0 1 2;
%      0 0 3 4 0 0 0 0;
%      5 6 0 0 0 0 0 0;
%      0 0 0 0 7 8 0 0;
%      0 0 9 10 0 0 0 0;
%      0 0 0 0 0 0 11 12;
%      0 0 0 0 13 14 0 0;
%      15 16 0 0 0 0 0 0];
% R = (1:1:16);
% [~,pos]=ismember(R,S);
% %yprob
% for ii = 1 : testBits-1
%     for jj = 1 : N
%         switch zkVec(jj,ii)
%             case 1+1i
%                 yprobTmp(:,:,(ii-1).*N+jj) = qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma).*qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma);
%             case 1-1i
%                 yprobTmp(:,:,(ii-1).*N+jj) = qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma).*( 1-qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma) );
%             case -1-1i
%                 yprobTmp(:,:,(ii-1).*N+jj) = ( 1-qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma) ).*( 1-qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma) );
%             case -1+1i
%                 yprobTmp(:,:,(ii-1).*N+jj) = ( 1-qfunc(-real(trellisOutput(8*jj-7 : 8*jj,:))./sigma) ).*qfunc(-imag(trellisOutput(8*jj-7 : 8*jj,:))./sigma);
%         end
% 
%         yprob(:,:,ii) = yprobTmp(:,:,(ii-1).*N+jj).*yprob(:,:,ii);
%     end
%     tmp = yprob(:,:,ii);
%     labelyj(ii,:) = tmp(pos);
%      labelyj(ii,:) =  labelyj(ii,:)./sum(labelyj(ii,:));
% end
% %%--------------------%
for idx = 1 : epoch
    for jj = 1 : iternum
        sumdLdRev = zeros(M,N);
        sumdLdImv = zeros(M,N);
        for ii = 1 : batch_size
            % [dLdRev(:,:,ii),dLdImv(:,:,ii),loss(ii),loss3(ii)] = GenGradient_SIMO(ResigVec( ii + (jj-1)*batch_size , : ),ImsigVec( ii + (jj-1)*batch_size , : ),Re_vMtr,Im_vMtr,labelVec( ii + (jj-1)*batch_size), labelyj( ii + (jj-1)*batch_size,:).');
            [dLdRev(:,:,ii),dLdImv(:,:,ii),loss(ii)] = GenGradient_SIMO(ResigVec( ii + (jj-1)*batch_size , : ),ImsigVec( ii + (jj-1)*batch_size , : ),Re_vMtr,Im_vMtr,labelVec( ii + (jj-1)*batch_size));            
            sumdLdRev = sumdLdRev + dLdRev(:,:,ii);
            sumdLdImv = sumdLdImv + dLdImv(:,:,ii);
        end
        gradRe = sumdLdRev ./ batch_size;
        gradIm = sumdLdImv ./ batch_size;
%         gradRe = mean(dLdRev,2);
%         gradIm = mean(dLdImv,2);
        lossAvg(jj) = mean(loss);
        % lossAvg2(jj) = mean(loss3);
        % Gradient Descent with Momentum
        vtRe = beta.* vtRe + (1 - beta).*(gradRe);
        vtIm = beta.* vtIm + (1 - beta).*(gradIm);

%         if jj>10 && ( abs(lossAvg(jj)-lossAvg(jj-10)) / abs(lossAvg(jj)) < 1e-5 )
%             break;
%         end
%         Re_vMtr = Re_vMtr - lr.*gradRe;
%         Im_vMtr = Im_vMtr - lr.*gradIm;
        Re_vMtr = Re_vMtr - lr.*vtRe;
        Im_vMtr = Im_vMtr - lr.*vtIm;
    end
    loss2(idx) = mean(lossAvg);
    % loss4(idx) = mean(lossAvg2);
end
% figure(2)
% subplot(1,2,2)
% plot(loss2,'ro-',LineWidth=1)
% xlabel('epoch')
% % xticks(0:4:101)
% ylabel('Cross Entropy')
% title('2 antennas')
% grid on
% hold on
% plot(loss4,'b+-',LineWidth=1)
% legend('QNN with Perfect CSI','QNN')
end