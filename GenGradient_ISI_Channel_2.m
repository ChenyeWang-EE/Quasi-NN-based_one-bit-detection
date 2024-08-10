function [dLdReh,dLdImh,loss2] = GenGradient_ISI_Channel_2(Resig,Imsig,Re_hVec,Im_hVec,label, xkMtr)

channelH_Vec = Re_hVec + Im_hVec*1i; % h/sigma
v = zeros(16,1);
for m = 1 : 16
    v(m) = channelH_Vec(1).*xkMtr(1,m) + channelH_Vec(2).*xkMtr(2,m);
end
Re_v = real(v);
Im_v = imag(v);

%% dxdReQ improved
% QRev_2 = qfunc(sign(Resig).*Re_v);
% QImv_2 = qfunc(sign(Imsig).*Im_v);
QRev = qfunc(-sign(Resig).*Re_v);
QImv = qfunc(-sign(Imsig).*Im_v);

% beta = log(QRev) + log(QImv);
% da_db = diag( exp(beta) ); % λ×λ

da_dRev = diag(  1/sqrt(2*pi) .* exp( log(QImv) ) .*exp(-Re_v.^2/2) .* sign(Resig)  ); % λ×λ
da_dImv = diag(  1/sqrt(2*pi) .* exp( log(QRev) ) .*exp(-Im_v.^2/2) .* sign(Imsig)  ); % λ×λ

% da_dRev = da_db * db_dRev;
% da_dImv = da_db * db_dImv; % λ×λ
dRev_dReh = zeros(16,2);
dRev_dImh = zeros(16,2);
dImv_dReh = zeros(16,2);
dImv_dImh = zeros(16,2);
for ii = 1 : 2
    % dRev_d(hR/simga)
    dRev_dReh(:,ii) = real(xkMtr(ii,:)).'; % λ×L
    % dRev_d(hI/simga)
    dRev_dImh(:,ii) = -imag(xkMtr(ii,:)).'; % λ×L
    % dImv_d(hR/simga)
    dImv_dReh(:,ii) = imag(xkMtr(ii,:)).'; % λ×L
    % dImv_d(hI/simga)
    dImv_dImh(:,ii) = real(xkMtr(ii,:)).'; % λ×L
end

da_dReh = da_dRev * dRev_dReh + da_dImv * dImv_dReh; % λ×L

da_dImh = da_dRev * dRev_dImh + da_dImv * dImv_dImh; % λ×L
%% compute dL_da (1×L)
sigma_Rev = QRev;
sigma_Imv = QImv;

a = exp(log(sigma_Rev) + log(sigma_Imv)); % λ×1
suma = sum(a);
prob_j = zeros(size(a));
for ii = 1 : length(a)
    prob_j(ii) =  a(ii) / suma ;
end

yj = one_hot_encode_ISI(label);
% yj = one_hot_encode(label); 
dLda = 1/suma.*( ones(size(yj)) - yj./prob_j ).'; % 1×λ

dLdReh = dLda * da_dReh; % 1×L
dLdImh = dLda * da_dImh; % 1×L

% da_d
% dxdReQ = QImv;
% dxdImQ = QRev;
% dReQ_dRev = -1/sqrt(2*pi).*exp(-Re_v.^2/2);
% % dReQ_dRev_2 = -sqrt(2).*Resig./sqrt(2*pi).*exp(-Re_v.^2/2);
% dx_dRev = dxdReQ .* dReQ_dRev;

% sigma_Rev = QRev;
% sigma_Imv = QImv;
% 
% 
% 
% dImQ_dImv = -1/sqrt(2*pi).*exp(-Im_v.^2/2) ;
% % dImQ_dImv_2 = -sqrt(2).*Imsig./sqrt(2*pi).*exp(-Im_v.^2/2) ;
% dx_dImv = dxdImQ .* dImQ_dImv;
% x = sigma_Rev.*sigma_Imv;
% sumx = sum(x);
% 
% prob_j = zeros(size(x));
% for ii = 1 : length(x)
%     prob_j(ii) =  x(ii) / sumx ;
% end
% 
% yj = one_hot_encode_ISI(label);
% % yj = one_hot_encode(label); 
% dLdx = 1/sumx.*( ones(size(yj)) - yj./prob_j );
% dLdRev = dLdx .* dx_dRev;
% dLdImv = dLdx .* dx_dImv;
% 
% dRev_dReh = zeros(16,2);
% dRev_dImh = zeros(16,2);
% dImv_dReh = zeros(16,2);
% dImv_dImh = zeros(16,2);
% dLdReh = zeros(1,2);
% dLdImh = zeros(1,2);
% for ii = 1 : 2
%     % dLd(hR/simga)
%     dRev_dReh(:,ii) = sign(Resig).*real(xkMtr(ii,:)).';
%     dRev_dImh(:,ii) = -sign(Resig).*imag(xkMtr(ii,:)).';
%     % dLd(hI/simga)
%     dImv_dReh(:,ii) = sign(Imsig).*imag(xkMtr(ii,:)).';
%     dImv_dImh(:,ii) = sign(Imsig).*real(xkMtr(ii,:)).';
% end
% 
% for ii = 1 : 2
%     dLdReh(ii) = sum( dLdRev.*dRev_dReh(:,ii) ) + sum( dLdImv.*dImv_dReh(:,ii) );
%     dLdImh(ii) = sum( dLdRev.*dRev_dImh(:,ii) ) + sum( dLdImv.*dImv_dImh(:,ii) );
%     if isnan(dLdReh(ii))
%         disp('NaN')
%         keyboard
%     end
% end
loss2 = -sum(yj.*log(prob_j));
loss2 = loss2 ./ length(yj);