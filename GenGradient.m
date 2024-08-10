function [dLdRev,dLdImv,loss] = GenGradient(Resig,Imsig,Re_v,Im_v,label,~)

QRev = qfunc(Re_v);
QImv = qfunc(Im_v);

dxdReQ =( ( 1/2+1/2*sign(Resig) ).*QRev.^(1/2*sign(Resig)-1/2).*(ones(size(QRev))-QRev).^(1/2-1/2*sign(Resig))...
    - (1/2-1/2*sign(Resig)).*(1-QRev).^(-1/2-1/2*sign(Resig)) ) .* (QImv.^(1/2+1/2*sign(Imsig))) .*(ones(size(QImv))-QImv).^(1/2-1/2*sign(Imsig));
dReQ_dRev = -1/sqrt(2*pi).*exp(-Re_v.^2/2) ;
dx_dRev = dxdReQ .* dReQ_dRev;

sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_v)) - QRev).^(1/2-1/2*sign(Resig));
sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_v)) - QImv).^(1/2-1/2*sign(Imsig));

dxdImQ =( ( 1/2+1/2*sign(Imsig) ).*QImv.^(1/2*sign(Imsig)-1/2).*(ones(size(QImv))-QImv).^(1/2-1/2*sign(Imsig))...
    - (1/2-1/2*sign(Imsig)).*(1-QImv).^(-1/2-1/2*sign(Imsig)) ) .* (QRev.^(1/2+1/2*sign(Resig))) .*(ones(size(QRev))-QRev).^(1/2-1/2*sign(Resig));
dImQ_dImv = -1/sqrt(2*pi).*exp(-Im_v.^2/2) ;
dx_dImv = dxdImQ .* dImQ_dImv;

x = sigma_Rev.*sigma_Imv;
sumx = sum(x);
prob_j = zeros(size(x));
for ii = 1 : length(x)
    prob_j(ii) =  x(ii) / sumx ;
end

yj = one_hot_encode_ISI(label);
% yj = one_hot_encode(label); 

dLdx = 1/sumx.*( ones(size(yj)) - yj./prob_j );
dLdRev = dLdx .* dx_dRev;
dLdImv = dLdx .* dx_dImv;
loss = -sum(yj.*log(prob_j));
loss = loss ./ length(yj);
% loss2 = -sum(yjtmp.*log(prob_j))./length(yjtmp);