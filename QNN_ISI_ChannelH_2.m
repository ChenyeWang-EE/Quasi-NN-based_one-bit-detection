function [prob_j] = QNN_ISI_ChannelH_2(testdata, Re_hVec,Im_hVec, xkMtr)

channelH_Vec = Re_hVec + Im_hVec*1i; % h/sigma
% v = -channelH .* xkVec;
v = zeros(16,1);
for m = 1 : 16
    v(m) = channelH_Vec(1).*xkMtr(1,m) + channelH_Vec(2).*xkMtr(2,m);
end
Re_v = real(v);
Im_v = imag(v);
Resig = real(testdata);
Imsig = imag(testdata);
QRev = qfunc(-sign(Resig).*Re_v);
QImv = qfunc(-sign(Imsig).*Im_v);
% sigma_Rev = QRev.^(1/2+1/2*sign(Resig)).*(ones(size(Re_v)) - QRev).^(1/2-1/2*sign(Resig));
% sigma_Imv = QImv.^(1/2+1/2*sign(Imsig)).*(ones(size(Im_v)) - QImv).^(1/2-1/2*sign(Imsig));
% x = sigma_Rev.*sigma_Imv;
% sumx = sum(x);
% prob_j = zeros(size(x));
% prob_j = prob_j.';
% for ii = 1 : length(x)
%     prob_j(ii) =  x(ii) / sumx ;
% end

% % QRev = qfunc(-sign(Resig).*Re_v);
% % QImv = qfunc(-sign(Imsig).*Im_v);
% % sigma_Rev = QRev;
% % sigma_Imv = QImv;
% % x = sigma_Rev.*sigma_Imv;
% % sumx = sum(x);
% % prob_j = zeros(size(x));
% % prob_j = prob_j.';
% % for ii = 1 : length(x)
% %     prob_j(ii) =  x(ii) / sumx ;
% % end

%% compute dL_da (1×L)
sigma_Rev = QRev;
sigma_Imv = QImv;

a = exp(log(sigma_Rev) + log(sigma_Imv)); % λ×1
suma = sum(a);
prob_j = zeros(size(a));
prob_j = prob_j.';
for ii = 1 : length(a)
    prob_j(ii) =  a(ii) / suma ;
end

end

