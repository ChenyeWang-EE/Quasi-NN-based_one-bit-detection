function [dLdReh2,dLdImh2,loss2] = GenGradient_SIMO_ISI_Channel(ResigVec,ImsigVec,Re_hMtr,Im_hMtr,label,xkMtr)

[N,~] = size(Re_hMtr);
M = 16;
channelH_Mtr = Re_hMtr + Im_hMtr.*1i;
vMtr = zeros(M,N);
for n = 1 : N
    for m = 1 : 16
        vMtr(m,n) = -channelH_Mtr(n,1).*xkMtr(1,m)-channelH_Mtr(n,2).*xkMtr(2,m);
    end
end
Re_vMtr = real(vMtr);
Im_vMtr = imag(vMtr);

%% v2
QRev2 = zeros(size(Re_vMtr));
QImv2 = zeros(size(Im_vMtr));
sigma_Rev2 = zeros(size(Re_vMtr));
sigma_Imv2 = zeros(size(Im_vMtr));
xn2 = zeros(size(Re_vMtr));
dxdReQ2 = zeros(size(Re_vMtr));
dReQ_dRev2 = zeros(size(Re_vMtr));
% dx_dRev = zeros(size(Re_vMtr));
dxdImQ2 = zeros(size(Im_vMtr));
dImQ_dImv2 = zeros(size(Im_vMtr));
% dx_dImv = zeros(size(Im_vMtr));
dLdRev2 = zeros(size(dReQ_dRev2));
dLdImv2 = zeros(size(dImQ_dImv2));
for n = 1 : N
    % QRev2(:,n) = qfunc(sign(ResigVec(n)).*Re_vMtr(:,n));
    QRev2(:,n) = qfunc( sqrt(2).*ResigVec(n).*Re_vMtr(:,n));
    % sqrt(2).*Resig.*
    % QImv2(:,n) = qfunc(sign(ImsigVec(n)).*Im_vMtr(:,n));
    QImv2(:,n) = qfunc(sqrt(2).*ImsigVec(n).*Im_vMtr(:,n));
    sigma_Rev2(:,n) = QRev2(:,n);
    sigma_Imv2(:,n) = QImv2(:,n);
    xn2(:,n) = sigma_Rev2(:,n).*sigma_Imv2(:,n);
end
x2 = prod(xn2,2);
for n = 1 : N
    dxdReQ2(:,n) = QImv2(:,n);
    dxdReQ2(:,n) = dxdReQ2(:,n) .* x2./xn2(:,n);
    dReQ_dRev2(:,n) = -1/sqrt(2*pi).*exp(-Re_vMtr(:,n).^2/2) ;

    dxdImQ2(:,n) = QRev2(:,n);
    dxdImQ2(:,n) = dxdImQ2(:,n) .* x2./xn2(:,n);
    dImQ_dImv2(:,n) = -1/sqrt(2*pi).*exp(-Im_vMtr(:,n).^2/2) ;
end
dx_dRev2 = dxdReQ2 .* dReQ_dRev2;
dx_dImv2 = dxdImQ2 .* dImQ_dImv2;
% sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_v)) - QRev).^(1/2-1/2*sign(Resig));
% sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_v)) - QImv).^(1/2-1/2*sign(Imsig));

% x = sigma_Rev.*sigma_Imv;
sumx2 = sum(x2);
prob_j2 = zeros(size(x2));
for ii = 1 : length(x2)
    prob_j2(ii) =  x2(ii) / sumx2 ;
end
yj = one_hot_encode_ISI(label);
% yj = one_hot_encode(label); 

dLdx2 = 1/sumx2.*( ones(size(yj)) - yj./prob_j2 );
for n = 1 : N
    dLdRev2(:,n) = dLdx2 .* dx_dRev2(:,n);
    dLdImv2(:,n) = dLdx2 .* dx_dImv2(:,n);
end

dRev_dReh2 = zeros(16,2,N);
dRev_dImh2 = zeros(16,2,N);
dImv_dReh2 = zeros(16,2,N);
dImv_dImh2 = zeros(16,2,N);
dLdReh2 = zeros(N,2);
dLdImh2 = zeros(N,2);
for n = 1 : N
    for ii = 1 : 2
        % dLd(hR/simga)
        dRev_dReh2(:,ii,n) = -sign(ResigVec(n)).*real(xkMtr(ii,:)).';
        dRev_dImh2(:,ii,n) = sign(ResigVec(n)).*imag(xkMtr(ii,:)).';
        % dLd(hI/simga)
        dImv_dReh2(:,ii,n) = -sign(ImsigVec(n)).*imag(xkMtr(ii,:)).';
        dImv_dImh2(:,ii,n) = -sign(ImsigVec(n)).*real(xkMtr(ii,:)).';
    end
end

for n = 1 : N
    for ii = 1 : 2
        dLdReh2(n,ii) = sum( dLdRev2(:,n).*dRev_dReh2(:,ii,n) ) + sum( dLdImv2(:,n).*dImv_dReh2(:,ii,n) );
        dLdImh2(n,ii) = sum( dLdRev2(:,n).*dRev_dImh2(:,ii,n) ) + sum( dLdImv2(:,n).*dImv_dImh2(:,ii,n) );
    end
end

loss2 = -sum(yj.*log(prob_j2))./length(yj);
% loss2 = -sum(yjtmp.*log(prob_j))./length(yjtmp);