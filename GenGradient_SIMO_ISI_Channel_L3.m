function [dLdReh,dLdImh,loss2] = GenGradient_SIMO_ISI_Channel_L3(ResigVec,ImsigVec,Re_hMtr,Im_hMtr,label,xkMtr)

[N,~] = size(Re_hMtr);
M = 64;
channelH_Mtr = Re_hMtr + Im_hMtr.*1i;
vMtr = zeros(M,N);
for n = 1 : N
    for m = 1 : M
        vMtr(m,n) = channelH_Mtr(n,1).*xkMtr(1,m) + channelH_Mtr(n,2).*xkMtr(2,m) + channelH_Mtr(n,3).*xkMtr(3,m);
    end
end
Re_vMtr = real(vMtr); % λ×M
Im_vMtr = imag(vMtr);

%% v2
QRev = zeros(size(Re_vMtr));
QImv = zeros(size(Im_vMtr));
sigma_Rev = zeros(size(Re_vMtr));
sigma_Imv = zeros(size(Im_vMtr));
% xn = zeros(size(Re_vMtr));
bn = zeros(size(Re_vMtr));
% dxdReQ = zeros(size(Re_vMtr));
% dReQ_dRev = zeros(size(Re_vMtr));
% dx_dRev = zeros(size(Re_vMtr));
% dxdImQ = zeros(size(Im_vMtr));
% dImQ_dImv = zeros(size(Im_vMtr));
% dx_dImv = zeros(size(Im_vMtr));
% dLdRev = zeros(size(dReQ_dRev));
% dLdImv = zeros(size(dImQ_dImv));


for n = 1 : N
    % QRev2(:,n) = qfunc(sign(ResigVec(n)).*Re_vMtr(:,n));
    QRev(:,n) = qfunc( -sign(ResigVec(n)).*Re_vMtr(:,n) );
    % sqrt(2).*Resig.*
    % QImv2(:,n) = qfunc(sign(ImsigVec(n)).*Im_vMtr(:,n));
    QImv(:,n) = qfunc( -sign(ImsigVec(n)).*Im_vMtr(:,n) );
    sigma_Rev(:,n) = QRev(:,n);
    sigma_Imv(:,n) = QImv(:,n);
    % xn(:,n) = sigma_Rev(:,n).*sigma_Imv(:,n);
    bn(:,n) = log( sigma_Rev(:,n) ) + log( sigma_Imv(:,n) );
end
% x = prod(xn,2);
b = sum(bn,2); % λ×1
a = exp( b ); % λ×1

da_db = diag( exp(b) );

da_dRev = zeros(M,M,N);
da_dImv = zeros(M,M,N);
for n = 1 : N
    da_dRev(:,:,n) = diag(  1/sqrt(2*pi) .* exp( b - bn(:,n) ) .* QImv(:,n) .*exp(-Re_vMtr(:,n).^2/2) .* sign(ResigVec(n))  ); % λ×λ
    da_dImv(:,:,n) = diag(  1/sqrt(2*pi) .* exp( b - bn(:,n) ) .* QRev(:,n) .*exp(-Im_vMtr(:,n).^2/2) .* sign(ImsigVec(n))  ); % λ×λ
end

dRev_dReh = zeros(64,3,N);
dRev_dImh = zeros(64,3,N);
dImv_dReh = zeros(64,3,N);
dImv_dImh = zeros(64,3,N);
dLdReh = zeros(N,3);
dLdImh = zeros(N,3);
for n = 1 : N
    for ii = 1 : 3
        % dRev_d(hR/simga)
        dRev_dReh(:,ii,n) = real(xkMtr(ii,:)).'; % λ×L
        % dRev_d(hI/simga)
        dRev_dImh(:,ii,n) = -imag(xkMtr(ii,:)).'; % λ×L
        % dImv_d(hR/simga)
        dImv_dReh(:,ii,n) = imag(xkMtr(ii,:)).'; % λ×L
        % dImv_d(hI/simga)
        dImv_dImh(:,ii,n) = real(xkMtr(ii,:)).'; % λ×L
    end
end
da_dReh = zeros(64,3,N);
da_dImh = zeros(64,3,N);
for n = 1 : N
    da_dReh(:,:,n) = da_dRev(:,:,n) * dRev_dReh(:,:,n) + da_dImv(:,:,n) * dImv_dReh(:,:,n); % λ×L
    da_dImh(:,:,n) = da_dRev(:,:,n) * dRev_dImh(:,:,n) + da_dImv(:,:,n) * dImv_dImh(:,:,n); % λ×L
end
%% compute dL_da (1×L)
suma = sum(a);
prob_j = zeros(size(a));
for ii = 1 : length(a)
    prob_j(ii) =  a(ii) / suma ;
end
yj = one_hot_encode_ISI_L3(label);
% yj = one_hot_encode(label); 

dLda = 1/suma.*( ones(size(yj)) - yj./prob_j ).'; % 1×λ

for n = 1 : N
    dLdReh(n,:)  = dLda * da_dReh(:,:,n); % 1×L
    dLdImh(n,:) = dLda * da_dImh(:,:,n); % 1×L
end

loss2 = -sum(yj.*log(prob_j))./length(yj);