function [dLdReh2,dLdImh2,loss2] = GenGradient_SIMO_Channel(ResigVec,ImsigVec,Re_hVec,Im_hVec,label,xkVec)

N = length(Re_hVec);  %antenna number
M = 8; % state
channelH_Vec = Re_hVec + Im_hVec*1i;
vMtr = zeros(M,N);
for nn = 1 : N
    vMtr(:,nn) = -channelH_Vec(nn) .* xkVec;
end
Re_vMtr = real(vMtr);
Im_vMtr = imag(vMtr);

% %% v1
% % [M,N] = size(Re_vMtr); %antenna number
% QRev = zeros(size(Re_vMtr));
% QImv = zeros(size(Im_vMtr));
% sigma_Rev = zeros(size(Re_vMtr));
% sigma_Imv = zeros(size(Im_vMtr));
% xn = zeros(size(Re_vMtr));
% dxdReQ = zeros(size(Re_vMtr));
% dReQ_dRev = zeros(size(Re_vMtr));
% % dx_dRev = zeros(size(Re_vMtr));
% dxdImQ = zeros(size(Im_vMtr));
% dImQ_dImv = zeros(size(Im_vMtr));
% % dx_dImv = zeros(size(Im_vMtr));
% dLdRev = zeros(size(dReQ_dRev));
% dLdImv = zeros(size(dImQ_dImv));
% 
% dLdReh = zeros(N,1);
% dLdImh = zeros(N,1);
% for n = 1 : N
%     QRev(:,n) = qfunc(Re_vMtr(:,n));
%     QImv(:,n) = qfunc(Im_vMtr(:,n));
%     sigma_Rev(:,n) = QRev(:,n).^(1/2+1/2*sign(ResigVec(n))).*(ones(M,1) - QRev(:,n)).^(1/2-1/2*sign(ResigVec(n)));
%     sigma_Imv(:,n) = QImv(:,n).^(1/2+1/2*sign(ImsigVec(n))).*(ones(M,1) - QImv(:,n)).^(1/2-1/2*sign(ImsigVec(n)));
%     xn(:,n) = sigma_Rev(:,n).*sigma_Imv(:,n);
% end
% x = prod(xn,2);
% for n = 1 : N
%     dxdReQ(:,n) =( ( 1/2+1/2*sign(ResigVec(n)) ).*QRev(:,n).^(1/2*sign(ResigVec(n))-1/2).*(ones(M,1)-QRev(:,n)).^(1/2-1/2*sign(ResigVec(n)))...
%         - (1/2-1/2*sign(ResigVec(n))).*(1-QRev(:,n)).^(-1/2-1/2*sign(ResigVec(n))) ).* (QImv(:,n).^(1/2+1/2*sign(ImsigVec(n)))) .*(ones(M,1)-QImv(:,n)).^(1/2-1/2*sign(ImsigVec(n)));
%     dxdReQ(:,n) = dxdReQ(:,n) .* x./xn(:,n);
%     dReQ_dRev(:,n) = -1/sqrt(2*pi).*exp(-Re_vMtr(:,n).^2/2) ;
% 
%     dxdImQ(:,n) =( ( 1/2+1/2*sign(ImsigVec(n)) ).*QImv(:,n).^(1/2*sign(ImsigVec(n))-1/2).*(ones(M,1)-QImv(:,n)).^(1/2-1/2*sign(ImsigVec(n)))...
%         - (1/2-1/2*sign(ImsigVec(n))).*(1-QImv(:,n)).^(-1/2-1/2*sign(ImsigVec(n))) ) .* (QRev(:,n).^(1/2+1/2*sign(ResigVec(n)))) .*(ones(M,1)-QRev(:,n)).^(1/2-1/2*sign(ResigVec(n)));
%     dxdImQ(:,n) = dxdImQ(:,n) .* x./xn(:,n);
%     dImQ_dImv(:,n) = -1/sqrt(2*pi).*exp(-Im_vMtr(:,n).^2/2) ;
% end
% dx_dRev = dxdReQ .* dReQ_dRev;
% dx_dImv = dxdImQ .* dImQ_dImv;
% % sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_v)) - QRev).^(1/2-1/2*sign(Resig));
% % sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_v)) - QImv).^(1/2-1/2*sign(Imsig));
% 
% % x = sigma_Rev.*sigma_Imv;
% sumx = sum(x);
% prob_j = zeros(size(x));
% for ii = 1 : length(x)
%     prob_j(ii) =  x(ii) / sumx ;
% end
% % yj = one_hot_encode_ISI(label); 
% yj = one_hot_encode(label); 
% 
% dLdx = 1/sumx.*( ones(size(yj)) - yj./prob_j );
% for n = 1 : N
%     dLdRev(:,n) = dLdx .* dx_dRev(:,n);
%     dLdImv(:,n) = dLdx .* dx_dImv(:,n);
% end
% %% dL/dReh; dL/dImh
% %dLd(hR/simga)
% dRev_dReh = -real(xkVec);
% dRev_dImh = imag(xkVec);
% %dLd(hI/simga)
% dImv_dReh = -imag(xkVec);
% dImv_dImh = -real(xkVec);
% 
% for n = 1 : N
%     dLdReh(n) = sum( dLdRev(:,n).*dRev_dReh ) + sum( dLdImv(:,n).*dImv_dReh );
%     dLdImh(n) = sum( dLdRev(:,n).*dRev_dImh ) + sum( dLdImv(:,n).*dImv_dImh );
% end
% 
% loss = -sum(yj.*log(prob_j))./length(yj);
% % loss2 = -sum(yjtmp.*log(prob_j))./length(yjtmp);

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

dLdReh2 = zeros(N,1);
dLdImh2 = zeros(N,1);
for n = 1 : N
    QRev2(:,n) = qfunc(sign(ResigVec(n)).*Re_vMtr(:,n));
    QImv2(:,n) = qfunc(sign(ImsigVec(n)).*Im_vMtr(:,n));
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
% yj = one_hot_encode_ISI(label); 
yj = one_hot_encode(label); 

dLdx2 = 1/sumx2.*( ones(size(yj)) - yj./prob_j2 );
for n = 1 : N
    dLdRev2(:,n) = dLdx2 .* dx_dRev2(:,n);
    dLdImv2(:,n) = dLdx2 .* dx_dImv2(:,n);
end
%% dL/dReh; dL/dImh
dRev_dReh2 = zeros(M,N);
dRev_dImh2 = zeros(M,N);
dImv_dReh2 = zeros(M,N);
dImv_dImh2 = zeros(M,N);
for n = 1 : N
    %dLd(hR/simga)
    dRev_dReh2(:,n) = -sign(ResigVec(n)).*real(xkVec);
    dRev_dImh2(:,n) = sign(ResigVec(n)).*imag(xkVec);
    %dLd(hI/simga)
    dImv_dReh2(:,n) = -sign(ImsigVec(n)).*imag(xkVec);
    dImv_dImh2(:,n) = -sign(ImsigVec(n)).*real(xkVec);
end

for n = 1 : N
    dLdReh2(n) = sum( dLdRev2(:,n).*dRev_dReh2(:,n) ) + sum( dLdImv2(:,n).*dImv_dReh2(:,n) );
    dLdImh2(n) = sum( dLdRev2(:,n).*dRev_dImh2(:,n) ) + sum( dLdImv2(:,n).*dImv_dImh2(:,n) );
end

loss2 = -sum(yj.*log(prob_j2))./length(yj);

end