function [prob_j] = QNN_SIMO_ISI_Channel_2(testdata, Re_hMtr,Im_hMtr,xkMtr)


% [M,N] = size(Re_vMtr);
% M = 16;
M = size(xkMtr,2);
[N,~] = size(Re_hMtr);
channelH_Mtr = Re_hMtr + Im_hMtr.*1i;
vMtr = zeros(M,N);
for n = 1 : N
    for m = 1 : M
        vMtr(m,n) = channelH_Mtr(n,1).*xkMtr(1,m) + channelH_Mtr(n,2).*xkMtr(2,m);
    end
end
Re_vMtr = real(vMtr);
Im_vMtr = imag(vMtr);
ResigVec = real(testdata);
ImsigVec = imag(testdata);

%% v2
QRev = zeros(size(Re_vMtr));
QImv = zeros(size(Im_vMtr));
sigma_Rev = zeros(size(Re_vMtr));
sigma_Imv = zeros(size(Im_vMtr));
% xn = zeros(size(sigma_Rev));
bn = zeros(size(Re_vMtr));
for n = 1 : N
    QRev(:,n) = qfunc( -sign(ResigVec(n)).*Re_vMtr(:,n));
    QImv(:,n) = qfunc( -sign(ImsigVec(n)).*Im_vMtr(:,n));
    sigma_Rev(:,n) = QRev(:,n);
    sigma_Imv(:,n) = QImv(:,n);
    % xn(:,n) = sigma_Rev(:,n).*sigma_Imv(:,n);
    bn(:,n) = log( sigma_Rev(:,n) ) + log( sigma_Imv(:,n) ); % λ×1
end
% x = prod(xn,2);
b = sum(bn,2); % λ×1
a = exp( b ); % λ×1
% sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_vMtr)) - QRev).^(1/2-1/2*sign(Resig));
% sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_vMtr)) - QImv).^(1/2-1/2*sign(Imsig));
% x = sigma_Rev.*sigma_Imv;
% sumx = sum(x);
% prob_j = zeros(size(x));
% for ii = 1 : length(x)
%     prob_j(ii) =  x(ii) / sumx ;
% end
suma = sum(a);
prob_j = zeros(size(a));
for ii = 1 : length(a)
    prob_j(ii) =  a(ii) / suma ;
end
end

