function [Re_h,Im_h,loss2] = forward_update_ISI_Channel_4QAM(ResigVec, ImsigVec, labelVec, beta, lr, batch_size,testBits)

% M = 16;
epoch = 100;
loss = zeros(1,batch_size);
% loss3 = zeros(1,batch_size);
Re_h = zeros(1,2);
Im_h = zeros(1,2);
% Re_v = -[0.7;0.7;-0.7;-0.7;0.7;-0.7;-0.7;0.7]/0.3162;
% Im_v = -[-0.7;0.7;0.7;-0.7;0.7;0.7;-0.7;-0.7]/0.3162;
iternum = floor((testBits-1) / batch_size);
lossAvg = zeros(1,iternum);
% lossAvg2 = zeros(1,iternum);
loss2 = zeros(1,epoch);
% loss4 = zeros(1,epoch);
vtRe = zeros(1,2); %initial gradient
vtIm = zeros(1,2);
% treliss output
map= [-1 / sqrt(2) + 1i / sqrt(2);
      -1 / sqrt(2) - 1i / sqrt(2);
       1 / sqrt(2) + 1i / sqrt(2);
       1 / sqrt(2) - 1i / sqrt(2)];
% xkMtr = [ exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(-0.25+0.5)), exp(1j*pi*(0.25+0.5)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)),...
%          exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25));
%          exp(1j*pi*(0-0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(0.5+0.25)),...
%          exp(1j*pi*(1-0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5+0.25)), exp(1j*pi*(-0.5+0.25))];
xkMtr = [ map(1), map(2), map(3), map(4), map(1), map(2), map(3), map(4), map(1), map(2), map(3), map(4), map(1), map(2), map(3), map(4);
          map(1), map(1), map(1), map(1), map(2), map(2), map(2), map(2), map(3), map(3), map(3), map(3), map(4), map(4), map(4), map(4)];
for idx = 1 : epoch
    for jj = 1 : iternum
        dLdReh = zeros(batch_size,2);
        dLdImh = zeros(batch_size,2);
        % dLdReh2 = zeros(batch_size,2);
        % dLdImh2 = zeros(batch_size,2);
        for ii = 1 : batch_size
            % [dLdRev(:,ii),dLdImv(:,ii),loss(ii),loss3(ii)] = GenGradient(ResigVec( ii + (jj-1)*batch_size ),ImsigVec( ii + (jj-1)*batch_size ),Re_v,Im_v,labelVec( ii + (jj-1)*batch_size),labelyj( ii + (jj-1)*batch_size,:).');
              
            % [dLdReh(ii,:),dLdImh(ii,:),loss(ii)] = GenGradient_ISI_Channel(ResigVec( ii + (jj-1)*batch_size ),ImsigVec( ii + (jj-1)*batch_size ),Re_h,Im_h,labelVec( ii + (jj-1)*batch_size), xkMtr);

              % [dLdReh(ii,:),dLdImh(ii,:),loss(ii)] = GenGradient_ISI_Channel_2(ResigVec( ii + (jj-1)*batch_size ),ImsigVec( ii + (jj-1)*batch_size ),Re_h,Im_h,labelVec( ii + (jj-1)*batch_size), xkMtr);
              [dLdReh(ii,:),dLdImh(ii,:),loss(ii)] = GenGradient_ISI_Channel_4QAM(ResigVec( ii + (jj-1)*batch_size ),ImsigVec( ii + (jj-1)*batch_size ),Re_h,Im_h,labelVec( ii + (jj-1)*batch_size), xkMtr);
              % if isnan(dLdReh(ii,1))
              %     disp('NaN')
              %     keyboard
              % end
        end
        gradRe = mean(dLdReh);
        gradIm = mean(dLdImh);
        lossAvg(jj) = mean(loss);
        % Gradient Descent with Momentum
        vtRe = beta.* vtRe + (1 - beta).*(gradRe); % beta is momentum term
        vtIm = beta.* vtIm + (1 - beta).*(gradIm);
%         if jj>10 && ( abs(lossAvg(jj)-lossAvg(jj-10)) / abs(lossAvg(jj-10)) < 1e-5 )
%             break;
%         end
        Re_h = Re_h - lr.*vtRe; % lr is learning rate
        Im_h = Im_h - lr.*vtIm;
    end
    loss2(idx) = mean(lossAvg);
end
% % plot
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