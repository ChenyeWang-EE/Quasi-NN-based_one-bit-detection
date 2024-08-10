function [Re_hMtr,Im_hMtr, loss2] = forward_update_SIMO_ISI_Channel(ResigVec, ImsigVec, labelVec, beta, lr, batch_size,testBits)

% M = 16;
[~,N] = size(ResigVec);
loss = zeros(1,batch_size);
% loss3 = zeros(1,batch_size);
% dLdRev = zeros(M,N,batch_size);
% dLdImv = zeros(M,N,batch_size);
% Re_vMtr = zeros(M,N);
% Im_vMtr = zeros(M,N);

dLdReh = zeros(N,2,batch_size);
dLdImh = zeros(N,2,batch_size);
Re_hMtr = zeros(N,2);
Im_hMtr = zeros(N,2);

% Re_v = -[0.7;0.7;-0.7;-0.7;0.7;-0.7;-0.7;0.7]/0.3162;
% Im_v = -[-0.7;0.7;0.7;-0.7;0.7;0.7;-0.7;-0.7]/0.3162;
iternum = floor((testBits-1) / batch_size);
lossAvg = zeros(1,iternum);
% lossAvg2 = zeros(1,iternum);
epoch = 100;
loss2 = zeros(1,epoch);
% loss4 = zeros(1,epoch);
vtRe = zeros(N,2); %initial gradient
vtIm = zeros(N,2);

xkMtr = [ exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(-0.25+0.5)), exp(1j*pi*(0.25+0.5)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)),...
         exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(-0.25-0.5)), exp(1j*pi*(0.25-0.5)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25));
         exp(1j*pi*(0-0.25)), exp(1j*pi*(0-0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0+0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5-0.25)), exp(1j*pi*(0.5+0.25)), exp(1j*pi*(0.5+0.25)),...
         exp(1j*pi*(1-0.25)), exp(1j*pi*(1-0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(1+0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5-0.25)), exp(1j*pi*(-0.5+0.25)), exp(1j*pi*(-0.5+0.25))];
for idx = 1 : epoch
    for jj = 1 : iternum
        sumdLdReh = zeros(N,2);
        sumdLdImh = zeros(N,2);
        for ii = 1 : batch_size
            % [dLdRev(:,:,ii),dLdImv(:,:,ii),loss(ii),loss3(ii)] = GenGradient_SIMO(ResigVec( ii + (jj-1)*batch_size , : ),ImsigVec( ii + (jj-1)*batch_size , : ),Re_vMtr,Im_vMtr,labelVec( ii + (jj-1)*batch_size), labelyj( ii + (jj-1)*batch_size,:).');
            % [dLdReh(:,:,ii),dLdImh(:,:,ii),loss(ii)] = GenGradient_SIMO_ISI_Channel(ResigVec( ii + (jj-1)*batch_size , : ),ImsigVec( ii + (jj-1)*batch_size , : ),Re_hMtr,Im_hMtr,labelVec( ii + (jj-1)*batch_size), xkMtr);            
            [dLdReh(:,:,ii),dLdImh(:,:,ii),loss(ii)] = GenGradient_SIMO_ISI_Channel_2(ResigVec( ii + (jj-1)*batch_size , : ),ImsigVec( ii + (jj-1)*batch_size , : ),Re_hMtr,Im_hMtr,labelVec( ii + (jj-1)*batch_size), xkMtr);            
            
            sumdLdReh = sumdLdReh + dLdReh(:,:,ii);
            sumdLdImh = sumdLdImh + dLdImh(:,:,ii);
        end
        gradRe = sumdLdReh ./ batch_size;
        gradIm = sumdLdImh ./ batch_size;
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
        Re_hMtr = Re_hMtr - lr.*vtRe;
        Im_hMtr = Im_hMtr - lr.*vtIm;
    end
    loss2(idx) = mean(lossAvg);
    % loss4(idx) = mean(lossAvg2);
end
% figure(2)
% subplot(1,2,2)
% plot(loss2,'r-',LineWidth=1)
% xlabel('epoch')
% % xticks(0:4:101)
% ylabel('Cross Entropy')
% title('2 antennas')
% grid on
% hold on
end