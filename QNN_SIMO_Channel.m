function [prob_j2] = QNN_SIMO_Channel(testdata, Re_hMtr,Im_hMtr,xkVec)


% [M,N] = size(Re_vMtr);
N = length(Re_hMtr);
M = 8;
channelH_Vec = Re_hMtr + 1i .* Im_hMtr;
vMtr = zeros(M,N);
for nn = 1 : N
    vMtr(:,nn) = -channelH_Vec(nn) .* xkVec;
end
Re_vMtr = real(vMtr);
Im_vMtr = imag(vMtr);

ResigVec = real(testdata);
ImsigVec = imag(testdata);
% %% v1
% QRev = zeros(size(Re_vMtr));
% QImv = zeros(size(Im_vMtr));
% sigma_Rev = zeros(size(Re_vMtr));
% sigma_Imv = zeros(size(Im_vMtr));
% xn = zeros(size(sigma_Rev));
% % QRev = qfunc(Re_vMtr);
% % QImv = qfunc(Im_vMtr);
% for n = 1 : N
%     QRev(:,n) = qfunc(Re_vMtr(:,n));
%     QImv(:,n) = qfunc(Im_vMtr(:,n));
%     sigma_Rev(:,n) = QRev(:,n).^(1/2+1/2*sign(ResigVec(n))).*(ones(M,1) - QRev(:,n)).^(1/2-1/2*sign(ResigVec(n)));
%     sigma_Imv(:,n) = QImv(:,n).^(1/2+1/2*sign(ImsigVec(n))).*(ones(M,1) - QImv(:,n)).^(1/2-1/2*sign(ImsigVec(n)));
%     xn(:,n) = sigma_Rev(:,n).*sigma_Imv(:,n);
% end
% x = prod(xn,2);
% 
% % sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_vMtr)) - QRev).^(1/2-1/2*sign(Resig));
% % sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_vMtr)) - QImv).^(1/2-1/2*sign(Imsig));
% % x = sigma_Rev.*sigma_Imv;
% sumx = sum(x);
% prob_j = zeros(size(x));
% for ii = 1 : length(x)
%     prob_j(ii) =  x(ii) / sumx ;
% end

%% v2
QRev2 = zeros(size(Re_vMtr));
QImv2 = zeros(size(Im_vMtr));
sigma_Rev2 = zeros(size(Re_vMtr));
sigma_Imv2 = zeros(size(Im_vMtr));
xn2 = zeros(size(sigma_Rev2));
% QRev = qfunc(Re_vMtr);
% QImv = qfunc(Im_vMtr);
for n = 1 : N
    QRev2(:,n) = qfunc(sign(ResigVec(n)).*Re_vMtr(:,n));
    QImv2(:,n) = qfunc(sign(ImsigVec(n)).*Im_vMtr(:,n));
    sigma_Rev2(:,n) = QRev2(:,n);
    sigma_Imv2(:,n) = QImv2(:,n);
    xn2(:,n) = sigma_Rev2(:,n).*sigma_Imv2(:,n);
end
x2 = prod(xn2,2);

% sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_vMtr)) - QRev).^(1/2-1/2*sign(Resig));
% sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_vMtr)) - QImv).^(1/2-1/2*sign(Imsig));
% x = sigma_Rev.*sigma_Imv;
sumx2 = sum(x2);
prob_j2 = zeros(size(x2));
for ii = 1 : length(x2)
    prob_j2(ii) =  x2(ii) / sumx2 ;
end
end

