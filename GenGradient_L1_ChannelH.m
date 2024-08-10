function [dLdReh2,dLdImh2,loss2] = GenGradient_L1_ChannelH(Resig,Imsig,Re_h,Im_h,label,xkVec)


channelH = Re_h + Im_h*1i;
v = -channelH .* xkVec;
Re_v = real(v);
Im_v = imag(v);

% %% v1
% QRev = qfunc(Re_v);
% QImv = qfunc(Im_v);
% 
% dxdReQ =( ( 1/2+1/2*sign(Resig) ).*QRev.^(1/2*sign(Resig)-1/2).*(ones(size(QRev))-QRev).^(1/2-1/2*sign(Resig))...
%     - (1/2-1/2*sign(Resig)).*(1-QRev).^(-1/2-1/2*sign(Resig)) ) .* (QImv.^(1/2+1/2*sign(Imsig))) .*(ones(size(QImv))-QImv).^(1/2-1/2*sign(Imsig));
% dReQ_dRev = -1/sqrt(2*pi).*exp(-Re_v.^2/2) ;
% dx_dRev = dxdReQ .* dReQ_dRev;
% 
% sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_v)) - QRev).^(1/2-1/2*sign(Resig));
% sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_v)) - QImv).^(1/2-1/2*sign(Imsig));
% 
% dxdImQ =( ( 1/2+1/2*sign(Imsig) ).*QImv.^(1/2*sign(Imsig)-1/2).*(ones(size(QImv))-QImv).^(1/2-1/2*sign(Imsig))...
%     - (1/2-1/2*sign(Imsig)).*(1-QImv).^(-1/2-1/2*sign(Imsig)) ) .* (QRev.^(1/2+1/2*sign(Resig))) .*(ones(size(QRev))-QRev).^(1/2-1/2*sign(Resig));
% dImQ_dImv = -1/sqrt(2*pi).*exp(-Im_v.^2/2) ;
% dx_dImv = dxdImQ .* dImQ_dImv;
% 
% x = sigma_Rev.*sigma_Imv;
% sumx = sum(x);
% prob_j = zeros(size(x));
% for ii = 1 : length(x)
%     prob_j(ii) =  x(ii) / sumx ;
% end
% 
% % yj = one_hot_encode_ISI(label); 
% yj = one_hot_encode(label); 
% 
% dLdx = 1/sumx.*( ones(size(yj)) - yj./prob_j );
% dLdRev = dLdx .* dx_dRev;
% dLdImv = dLdx .* dx_dImv;
% 
% %dLd(hR/simga)
% dRev_dReh = -real(xkVec);
% dRev_dImh = imag(xkVec);
% %dLd(hI/simga)
% dImv_dReh = -imag(xkVec);
% dImv_dImh = -real(xkVec);
% 
% dLdReh = sum( dLdRev.*dRev_dReh ) + sum( dLdImv.*dImv_dReh );
% dLdImh = sum( dLdRev.*dRev_dImh ) + sum( dLdImv.*dImv_dImh );
% 
% loss = -sum(yj.*log(prob_j));
% loss = loss ./ length(yj);
% % loss2 = -sum(yjtmp.*log(prob_j))./length(yjtmp);

%% v2
QRev2 = qfunc(sign(Resig).*Re_v);
QImv2 = qfunc(sign(Imsig).*Im_v);

dxdReQ2 = QImv2;
dReQ_dRev2 = -1/sqrt(2*pi).*exp(-Re_v.^2/2) ;
dx_dRev2 = dxdReQ2 .* dReQ_dRev2;

sigma_Rev2 = QRev2;
sigma_Imv2 = QImv2;

dxdImQ2 = QRev2;
dImQ_dImv2 = -1/sqrt(2*pi).*exp(-Im_v.^2/2) ;
dx_dImv2 = dxdImQ2 .* dImQ_dImv2;

x2 = sigma_Rev2.*sigma_Imv2;
sumx2 = sum(x2);
prob_j2 = zeros(size(x2));
for ii = 1 : length(x2)
    prob_j2(ii) =  x2(ii) / sumx2 ;
end

% yj = one_hot_encode_ISI(label); 
yj = one_hot_encode(label); 

dLdx2 = 1/sumx2.*( ones(size(yj)) - yj./prob_j2 );
dLdRev2 = dLdx2 .* dx_dRev2;
dLdImv2 = dLdx2 .* dx_dImv2;

%dLd(hR/simga)
dRev_dReh2 = -sign(Resig).*real(xkVec);
dRev_dImh2 = sign(Resig).*imag(xkVec);
%dLd(hI/simga)
dImv_dReh2 = -sign(Imsig).*imag(xkVec);
dImv_dImh2 = -sign(Imsig).*real(xkVec);

dLdReh2 = sum( dLdRev2.*dRev_dReh2 ) + sum( dLdImv2.*dImv_dReh2 );
dLdImh2 = sum( dLdRev2.*dRev_dImh2 ) + sum( dLdImv2.*dImv_dImh2 );

loss2 = -sum(yj.*log(prob_j2));
loss2 = loss2 ./ length(yj);

end