function [x2] = QNN_L1_ChannelH(testdata, Re_h,Im_h, xkVec)

channelH = Re_h + 1i*Im_h;
v = -channelH .* xkVec;
Re_v = real(v);
Im_v = imag(v);
Resig = real(testdata);
Imsig = imag(testdata);
% %% v1
% QRev = qfunc(Re_v);
% QImv = qfunc(Im_v);
% sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_v)) - QRev).^(1/2-1/2*sign(Resig));
% sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_v)) - QImv).^(1/2-1/2*sign(Imsig));
% x = sigma_Rev.*sigma_Imv;
% sumx = sum(x);
% prob_j = zeros(size(x));
% prob_j = prob_j.';
% for ii = 1 : length(x)
%     prob_j(ii) =  x(ii) / sumx ;
% end

%% v2
QRev2 = qfunc(sign(Resig).*Re_v);
QImv2 = qfunc(sign(Imsig).*Im_v);
sigma_Rev2 = QRev2;
sigma_Imv2 = QImv2;
x2 = sigma_Rev2.*sigma_Imv2;
% sumx2 = sum(x2);
% prob_j2 = zeros(size(x2));
% prob_j2 = prob_j2.';
% for ii = 1 : length(x2)
%     prob_j2(ii) =  x2(ii) / sumx2 ;
% end
end

