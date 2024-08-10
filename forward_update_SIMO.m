function [Re_vMtr,Im_vMtr, loss] = forward_update_SIMO(ResigVec, ImsigVec, labelVec, beta, lr, batch_size,testBits)

M = 8;
[~,N] = size(ResigVec);
loss = zeros(1,batch_size);
dLdRev = zeros(M,N,batch_size);
dLdImv = zeros(M,N,batch_size);
Re_vMtr = zeros(M,N);
Im_vMtr = zeros(M,N);
% Re_v = -[0.7;0.7;-0.7;-0.7;0.7;-0.7;-0.7;0.7]/0.3162;
% Im_v = -[-0.7;0.7;0.7;-0.7;0.7;0.7;-0.7;-0.7]/0.3162;
iternum = floor(testBits / batch_size);
lossAvg = zeros(1,iternum);
epoch = 50;
vtRe = zeros(M,N); %initial gradient
vtIm = zeros(M,N);
for idx = 1 : epoch
    for jj = 1 : iternum
        sumdLdRev = zeros(M,N);
        sumdLdImv = zeros(M,N);
        for ii = 1 : batch_size
            [dLdRev(:,:,ii),dLdImv(:,:,ii),loss(ii)] = GenGradient_SIMO(ResigVec( ii + (jj-1)*batch_size , : ),ImsigVec( ii + (jj-1)*batch_size , : ),Re_vMtr,Im_vMtr,labelVec( ii + (jj-1)*batch_size));
            sumdLdRev = sumdLdRev + dLdRev(:,:,ii);
            sumdLdImv = sumdLdImv + dLdImv(:,:,ii);
        end
        gradRe = sumdLdRev ./ batch_size;
        gradIm = sumdLdImv ./ batch_size;
%         gradRe = mean(dLdRev,2);
%         gradIm = mean(dLdImv,2);
        lossAvg(jj) = mean(loss);
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
end

end