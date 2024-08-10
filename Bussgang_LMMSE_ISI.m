function hHeadBLMMSE = Bussgang_LMMSE_ISI(xk,z,antennaNum)

x1 = xk(2:end);
x2 = xk(1:end-1);
x = [x1; x2];

[~,xLen] =size(x);
zvec = z(:);
s = kron(x.',eye(antennaNum));
Cy = kron(x.'*conj(x),eye(antennaNum))+eye(antennaNum*xLen);
ss = diag( Cy ).^(-1/2) ;
A = sqrt(2/pi).*real(diag(ss));
ATilde = A*s;
% BLMMSE
Cz = 2/pi.*( asin( real(diag(ss))*real(Cy)*real(diag(ss)) ) + 1i*asin( real(diag(ss))*imag(Cy)*real(diag(ss)) ) );
hHeadBLMMSEtmp = ATilde'*(Cz\zvec);
hHeadBLMMSE = reshape(hHeadBLMMSEtmp,[antennaNum,length(hHeadBLMMSEtmp)/antennaNum]);
% BLS
hHeadBLStmp =  (ATilde'*ATilde)\ATilde'*zvec;
hHeadBLS = reshape(hHeadBLStmp,[antennaNum,length(hHeadBLStmp)/antennaNum]);
end