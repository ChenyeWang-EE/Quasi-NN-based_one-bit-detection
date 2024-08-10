function [dLdReh,dLdImh,loss2] = GenGradient_ISI_Channel_4QAM_L3(Resig,Imsig,Re_hVec,Im_hVec,label, xkMtr)

channelH_Vec = Re_hVec + Im_hVec*1i; % h/sigma
v = zeros(64,1);
for m = 1 : 64
    v(m) = channelH_Vec(1).*xkMtr(1,m) + channelH_Vec(2).*xkMtr(2,m) + channelH_Vec(3).*xkMtr(3,m);
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
dRev_dReh = zeros(64,3);
dRev_dImh = zeros(64,3);
dImv_dReh = zeros(64,3);
dImv_dImh = zeros(64,3);
for ii = 1 : 3
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

yj = one_hot_encode_ISI_L3(label);
% yj = one_hot_encode(label); 
dLda = 1/suma.*( ones(size(yj)) - yj./prob_j ).'; % 1×λ

dLdReh = dLda * da_dReh; % 1×L
dLdImh = dLda * da_dImh; % 1×L


loss2 = -sum(yj.*log(prob_j));
loss2 = loss2 ./ length(yj);