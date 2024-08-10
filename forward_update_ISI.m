function [Re_v,Im_v,loss] = forward_update_ISI(ResigVec, ImsigVec, labelVec, beta, lr, batch_size,testBits)

M = 16;
epoch = 50;
loss = zeros(1,batch_size);
% loss3 = zeros(1,batch_size);
dLdRev = zeros(M,batch_size);
dLdImv = zeros(M,batch_size);
Re_v = zeros(M,1);
Im_v = zeros(M,1);
% Re_v = -[0.7;0.7;-0.7;-0.7;0.7;-0.7;-0.7;0.7]/0.3162;
% Im_v = -[-0.7;0.7;0.7;-0.7;0.7;0.7;-0.7;-0.7]/0.3162;
iternum = floor((testBits-1) / batch_size);
lossAvg = zeros(1,iternum);
% lossAvg2 = zeros(1,iternum);
loss2 = zeros(1,epoch);
% loss4 = zeros(1,epoch);
vtRe = zeros(M,1); %initial gradient
vtIm = zeros(M,1);
% %%---------------------%%
% labelyj = zeros(testBits-1,16);
% yprob = zeros(8,8,testBits-1);
% trellisOutput = statechangeoutput(channelH);
% %yprob
% zk = ResigVec + 1i*ImsigVec;
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
% for ii = 1 : testBits-1
%     switch zk(ii)
%         case 1+1i
%             yprob(:,:,ii) = qfunc(-real(trellisOutput)./sigma).*qfunc(-imag(trellisOutput)./sigma);
%         case 1-1i
%             yprob(:,:,ii) = qfunc(-real(trellisOutput)./sigma).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1-1i
%             yprob(:,:,ii) = ( 1-qfunc(-real(trellisOutput)./sigma) ).*( 1-qfunc(-imag(trellisOutput)./sigma) );
%         case -1+1i
%             yprob(:,:,ii) = ( 1-qfunc(-real(trellisOutput)./sigma) ).*qfunc(-imag(trellisOutput)./sigma);
%     end
%     tmp = yprob(:,:,ii);
%     labelyj(ii,:) = tmp(pos);
%     labelyj(ii,:) = labelyj(ii,:)./ sum(labelyj(ii,:));
% end
% %%--------------------%

for idx = 1 : epoch
    for jj = 1 : iternum
        for ii = 1 : batch_size
            % [dLdRev(:,ii),dLdImv(:,ii),loss(ii),loss3(ii)] = GenGradient(ResigVec( ii + (jj-1)*batch_size ),ImsigVec( ii + (jj-1)*batch_size ),Re_v,Im_v,labelVec( ii + (jj-1)*batch_size),labelyj( ii + (jj-1)*batch_size,:).');
              [dLdRev(:,ii),dLdImv(:,ii),loss(ii)] = GenGradient(ResigVec( ii + (jj-1)*batch_size ),ImsigVec( ii + (jj-1)*batch_size ),Re_v,Im_v,labelVec( ii + (jj-1)*batch_size));
        end
        gradRe = mean(dLdRev,2);
        gradIm = mean(dLdImv,2);
        lossAvg(jj) = mean(loss);
        % lossAvg2(jj) = mean(loss3);
        % Gradient Descent with Momentum
        vtRe = beta.* vtRe + (1 - beta).*(gradRe);
        vtIm = beta.* vtIm + (1 - beta).*(gradIm);
%         if jj>10 && ( abs(lossAvg(jj)-lossAvg(jj-10)) / abs(lossAvg(jj-10)) < 1e-5 )
%             break;
%         end
        Re_v = Re_v - lr.*vtRe;
        Im_v = Im_v - lr.*vtIm;
    end
    loss2(idx) = mean(lossAvg);
    % loss4(idx) = mean(lossAvg2);
end
%plot
% figure(2)
% subplot(1,2,1)
% plot(loss2,'ro-',LineWidth=1)
% xlabel('epoch')
% % xticks(0:4:101)
% ylabel('Cross Entropy')
% title('Single antenna')
% grid on
% hold on 
% plot(loss4,'b+-',LineWidth=1)
% legend('QNN with Perfect CSI','QNN')
% sgtitle('SNR = 100dB')
end